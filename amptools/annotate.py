import sys
import itertools
import logging; log = logging.getLogger(__name__)
import pysam
import json
from collections import Counter

import amplicon
import stats

TAG_COUNT = 'mc'
TAG_AMP = 'ea'
TAG_BC = 'BC'
TAG_RG = 'RG'


def _read_trim_file(trim_file):
    """ read barcodes or molecular counters from a cutadapt trim file
        returns a dictionary of (accession, sequence)
    """
    # the read part is the everything before the first space
    # the accession is everything afterwards
    log.info('reading file {0}'.format(trim_file))
    read_mids = itertools.imap(
        lambda line: line.rstrip().split(' ',1),
        file(trim_file)
    )

    try:
        read_mids = dict(((y, x) for (x,y) in read_mids))
    except Exception, e:
        raise Exception('could not parse file %s: %s' % (trim_file, e))

    return read_mids


class MidAnnotator(object):
    """ Annotate BAM file with read groups (RGs) based on molecular barcodes.

        This annotator adds RGs to the header and assigns each read a RG tag with the
        reag group and a BC tag with the barcode read.
    """

    # TODO: if the BAM is in the same order as the reads, can we avoid preloading read mids
    # TODO: write the barcode quality attribute

    @classmethod
    def customize_parser(cls, parser):
        group = parser.add_argument_group('RG annotation', cls.__doc__)
        group.add_argument('--rgs', type=str, help='file containing whitespace separated BC, RG pairs (enables annotator)')
        group.add_argument('--bcs-read', type=str, help='file containing whitespace separated BC, read accession pairs')
        group.add_argument('--rgs-read', type=str, help='file containing whitespace separated RG, read accession pairs')
        group.add_argument('--library', type=str, help='(optional) library to use in RG header')
        group.add_argument('--platform', type=str, help='(optional) platform to use in RG header')
        group.add_argument('--exclude-rg', type=str, help='(optional) platform to use in RG header')
        group.add_argument('--offbyone', action='store_true', help='Allow off by one errors in the MID')


    def __init__(self, args, header):

        self.counts = Counter()

        assert args.rgs, 'Need --rgs'
        assert args.bcs_read or args.rgs_read, 'please provide either --rgs-read or --bcs-read'

        self.exclude = args.exclude_rg

        # parse the trim file of actually read mids
        if args.bcs_read:
            self.read_bcs = _read_trim_file(args.bcs_read)
            self.read_rgs = None
        else:
            self.read_rgs = _read_trim_file(args.rgs_read)
            self.read_bcs = None

        mids =  itertools.imap(
            lambda line: line.rstrip().split('\t', 1),
            file(args.rgs)
        )

        self.mids = dict(mids)
        log.info('read {0} mids'.format(len(self.mids)))
        for x in self.mids.values():
            self.counts[x] = 0
        # update the header
        RGS = []
        template = {}
        if args.platform:
            template['PL'] = args.platform
        if args.library:
            template['LB'] = args.library
        for mid in self.mids.values():
            if not mid == self.exclude:
                entry = {'ID': mid, 'SM': mid}
                entry.update(template)
                RGS.append(entry)

        header['RG'] = RGS

        # precompute off by one errors:
        if args.offbyone:
            for bc, mid in self.mids.items():
                for (i, base) in enumerate(bc):
                    for sub in 'ACGT':
                        if sub != base:
                            alt = bc[:i] + sub + bc[i+1:]
                            if alt not in self.mids:
                                self.mids[alt] = bc


    def match_read(self, read_mid):
        return self.mids[read_mid]

    def __call__(self, read):
        try:
            if self.read_bcs:
                read_bc = self.read_bcs[read.qname]
                log.info('rbc ' + read_bc + ' ' + ','.join(self.mids.keys()))
                RG = self.match_read(read_bc)
                etags = [('RG', RG), ('BC', read_bc)]
            else:
                RG = self.read_rgs[read.qname]
                etags = [('RG', RG)]

        except KeyError:
            self.counts['No match'] += 1
            #TODO: fallback matching and stats
            return

        self.counts[RG] += 1
        read.tags = read.tags + etags
        if RG == self.exclude:
            return False
        return read

    def report(self):
        print 'sample reads'

        for k,v in sorted(self.counts.items()):
           print k, v
        pass


class DbrAnnotator(object):
    """ Annotate BAM file with molecular counters (MCs).

        This annotator adds a MC tag for each read contaning any molecular
        counter sequence read.
    """

    # TODO: same as MidAnnotator

    @classmethod
    def customize_parser(cls, parser):
        group = parser.add_argument_group('MC annotation', cls.__doc__)
        group.add_argument('--counters', type=str,
                help='File containing whitespace separated MC, read accesion')

    def __init__(self, args, header):
        self.counts = Counter()
        self.read_mids = _read_trim_file(args.counters)

    def __call__(self, read):
        try:
            MC = self.read_mids[read.qname]
            if MC != '':
                read.tags = read.tags + [(TAG_COUNT, MC)]
            self.counts[MC] += 1
            return read
        except KeyError:
            self.counts[None] += 1
            # TODO: stats for missing counter

    def report(self):
        print 'counter reads'

        for k,v in sorted(self.counts.items()):
           print k, v
        pass



class AmpliconAnnotator(object):
    """ Annotate reads that match expected amplicons (EAs).

        Mark each read that matches and expected amplicon with an EA tag and
        add a EA section to the header.
    """
    @classmethod
    def customize_parser(cls, parser):
        group = parser.add_argument_group('EA annotation', cls.__doc__)

        group.add_argument('--amps', type=str,
                help='Delimited file describing amplicons (required)')
        group.add_argument('--id-column', type=str,
                help='ID column (default id)', default='id')
        group.add_argument('--amplicon-column', type=str,
                help='Amplicon coordinate column (default amplicon)', default='amplicon')
        group.add_argument('--trim-column', type=str,
                help='Amplicon trim coordinates (default trim)', default='trim')
        group.add_argument('--delimiter', type=str,
                help='file delimiter (default TAB)', default='\t')
        group.add_argument('--offset-allowed', type=int,
                help='Allowed bases between read start and amplicon start (default 10)', default=10)
        group.add_argument('--exclude-offtarget', action='store_true', default=False,
                help='Exclude reads not matching target amplicons')
        # FIXME: remove argument?
        group.add_argument('--clip', action='store_true')


    def __init__(self, args, header):
        self.stats = stats.Stats('')
        self.amplicons = amplicon.load_amplicons(args.amps, self.stats, args)
        self.clip = args.clip
        self.exclude_offtarget = args.exclude_offtarget

        AMS = []
        for amp in self.amplicons:
            AMS.append(json.dumps({
                'type': 'ea',
                'id': amp.external_id,
                'ac': '%s:%s-%s' % (amp.chr, amp.start, amp.end),
                'tc': '%s:%s-%s' % (amp.chr, amp.trim_start, amp.trim_end),
                'st': str(amp.strand)
                }))

        header['CO'] = header.get('CO', []) + AMS

        # create a list of lists ref by tid
        self._amps_by_chr = []
        for _ in range(args.input.nreferences):
            self._amps_by_chr.append([])

        for a in self.amplicons:
            self._amps_by_chr[args.input.gettid(a.chr)].append(a)

    def __call__(self, read):
        for amp in self._amps_by_chr[read.tid]:
            if amp.matches(read):
                # FIXME: amplicon mark method
                read.tags = read.tags + [(TAG_AMP, amp.external_id)]
                if self.clip:
                    clipped = amp.clip(read)
                    if clipped:
                        return clipped
                    else:
                        return False
                else:
                    return read
        if self.exclude_offtarget:
            return False
        return read


    def report(self):
        self.stats.report(sys.stdout)


def annotate(args):
    """ Annotate reads in a SAM file with tags.

        Use one or more available annotators below to add tags to a SAM file.

    """
    inp = args.input = pysam.Samfile(args.input)

    header = inp.header
    annotators = []

    try:
        if args.amps:
            annotators.append(AmpliconAnnotator(args, header))
        if args.rgs or args.rgs_read:
            annotators.append(MidAnnotator(args, header))
        if args.counters:
            annotators.append(DbrAnnotator(args, header))
    except Exception, e:
        print e
        raise
        sys.exit(1)

    assert 'SQ' in header # http://code.google.com/p/pysam/issues/detail?id=84
    oup = pysam.Samfile(args.output, 'wb', header=header)

    log.info('begin read annotation')
    for read in inp:
        include = True
        for annotator in annotators:
            # Annotators return False to exclude
            if annotator(read) is False:
                include = False
                break
        if include:
            oup.write(read)


    for a in annotators:
        a.report()


def duplicates(args):
    """ Mark duplicates using a molecular counter.

        The file should contain MC tags.  Duplicates are detected by looking
        for reads in the same start position and orientation and with the same MC tag.

        WARNING: this unsorts your input
    """
    inp = pysam.Samfile(args.input, "rb" )
    outp = pysam.Samfile(args.output, "wb", template=inp)

    # TODO: check sorted
    # TODO: report duplicates/unannotated stats
    # TODO: add sort option when output set

    current_pos = None
    to_check = {}

    def filter_dups(group, positions):
        """ filter a set of start positions for the best read by mapping quality """
        entries = []
        for pos in current:
            log.debug('using entries from %s, %s' % (group, pos))
            entries.extend(to_check.pop((group[0], group[1], pos)))

        # hash based on the RG and the DBR
        log.debug('filter %s reads' % len(entries))
        groups = {}
        for e in entries:
            try:
                rg, dbr = e.opt('RG'), e.opt(TAG_COUNT)
                # FIXME: parameterize this DBR validation code
                #if len(dbr) != 2 or 'N' in dbr:
                #    continue

                group = (rg, dbr)
                groups[group] = groups.get(group, []) + [e]
            except KeyError:
                log.debug('read %s missing required tags' % e.qname)


        # return the best read in each group
        for dups in groups.values():

            # FIXME: option to choose best read stratedy
            keyfunc = lambda x: x.mapq
            #if opts.length:
            #    keyfunc = lambda x: x.qlen

            ordered = sorted(dups, key=keyfunc)
            for not_best in ordered[:-1]:
                not_best.is_duplicate = True
                yield not_best
            # mark duplicates

            yield ordered[-1]


    n=0
    for (i, entry) in enumerate(inp):
        if (i % 10) == 0:
            log.info('processed %(i)s reads' % locals())

        # TODO: downsample?
        #if args.random and random.random() > args.random:
        #    continue

        # find out the start and orientation of the read
        is_reverse = entry.is_reverse
        start = entry.pos if not is_reverse else entry.aend

        # append to the to_check index
        position_index = (entry.rname, is_reverse, start)
        to_check[position_index] = to_check.get(position_index, []) + [entry]

        if entry.pos != current_pos:
            current_pos = entry.pos

            # TODO: could process some of the reads when the reference changes
            # if entry.rname != rname: PASS

    merge_distance = 4
    indexes = sorted(to_check.keys())
    log.debug(indexes)
    for group, inds in itertools.groupby(indexes, key=lambda x: x[:2]):

        log.debug('processing group %(group)s' % locals())
        positions = [x[2] for x in inds]
        current = [positions.pop(0)]

        while True:
            head = current[-1]
            if len(positions):
                consider = positions[0]
                log.debug('head is %(head)s, considering %(consider)s' % locals())
                if abs(consider - head) <= merge_distance:
                    current.append(positions.pop(0))
                    continue

            log.debug('deduping %(current)s' % locals())
            for nondup in filter_dups(group, positions):
                outp.write(nondup)
            if len(positions):

                current = [positions.pop(0)]
            else:
                break

        #for nondup in filter_dups(group, current):
        #    log.debug('writing read')
        #    outp.write(nondup)
#



