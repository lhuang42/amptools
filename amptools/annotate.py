import sys
import itertools
import random
from collections import Counter, defaultdict

import pysam

import amplicon
import stats


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
        group.add_argument('--mids', type=str, help='file containing whitespace separated BC, RG pairs')
        group.add_argument('--mids-read', type=str, help='file containing whitespace separated BC, read accession pairs')
        group.add_argument('--library', type=str, help='(optional) library to use in RG header')
        group.add_argument('--platform', type=str, help='(optional) platform to use in RG header')


    def __init__(self, args, header):

        self.counts = Counter()

        assert args.mids and args.mids_read, 'please provide both --mids and --mids-read'

        # parse the trim file of actually read mids
        read_mids = itertools.imap(
            lambda line: line.rstrip().split('\t', 1),
            file(args.mids_read)
        )

        self.read_mids = dict(read_mids)

        # get the expected mids
        mids =  itertools.imap(
            lambda line: line.rstrip().split(),
            file(args.mids)
        )


        self.mids = dict(mids)
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
            entry = {'ID': mid, 'SM': mid}
            entry.update(template)
            RGS.append(entry)

        header['RG'] = RGS


    def match_read(self, read_mid):
        return self.mids[read_mid]

    def __call__(self, read):
        try:
            read_mid = self.read_mids[read.qname]
            RG = self.match_read(read_mid)
        except KeyError:
            self.counts['No match'] += 1
            #TODO: fallback matching and stats
            return

        self.counts[RG] += 1
        read.tags = read.tags + [('RG', RG), ('BC', read_mid)]

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
        # TODO: factor out reading cutadapt data
        # parse the trim file of actually read mids
        read_mids = itertools.imap(
            lambda line: line.rstrip().split('\t', 1),
            file(args.counters)
        )

        try:
            self.read_mids = dict(read_mids)
        except Exception, e:
            raise Exception('could not parse counter file %s: %s' % (args.counters, e))


    def __call__(self, read):
        try:
            MC = self.read_mids[read.qname]
            read.tags = read.tags + [('MC', MC)]
            self.counts[MC] += 1
        except KeyError:
            self.counts[None] += 1
            #print 'missing', read.qname
            pass
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
        # FIXME: remove argument?
        group.add_argument('--clip', action='store_true')


    def __init__(self, args, header):
        self.stats = stats.Stats('')
        self.amplicons = amplicon.load_amplicons(args.amps, self.stats, args)

        self.clip = args.clip

        AMS = []
        for amp in self.amplicons:
            AMS.append({
                'ID': amp.external_id,
                'AC': '%s:%s-%s' % (amp.chr, amp.start, amp.end),
                'TC': '%s:%s-%s' % (amp.chr, amp.trim_start, amp.trim_end),
                'ST': str(amp.strand)
                })

        header['EA'] = AMS

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
                read.tags = read.tags + [('EA', amp.external_id)]
                if self.clip:
                    amp.clip(read)
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
        if args.mids:
            annotators.append(MidAnnotator(args, header))
        if args.counters:
            annotators.append(DbrAnnotator(args, header))
    except Exception, e:
        print e
        sys.exit(1)

    assert 'SQ' in header # http://code.google.com/p/pysam/issues/detail?id=84
    oup = pysam.Samfile(args.output, 'wb', header=header)

    for read in inp:
        for annotator in annotators:
            annotator(read)
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
            entries.extend(to_check.pop((group[0], group[1], pos)))

        # hash based on the RG and the DBR
        groups = {}
        for e in entries:
            rg, dbr = e.opt('RG'), e.opt('MC')

            # FIXME: parameterize this DBR validation code
            if len(dbr) != 3 or 'N' in dbr:
                continue

            group = (rg, dbr)
            groups[group] = groups.get(group, []) + [e]

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
        if (i % 10000) == 0:
            print('processed %(i)s reads' % locals())

        #if args.random and random.random() > args.random:
        #    continue
        n+=1
        is_reverse = entry.is_reverse
        start = entry.pos if not is_reverse else entry.aend

        position_index = (entry.rname, is_reverse, start)
        to_check[position_index] = to_check.get(position_index, []) + [entry]

        if entry.pos != current_pos:
            current_pos = entry.pos

            # could process some of the reads when the reference changes
            # if entry.rname != rname: PASS

            # find positions near enough to be the same read and hand off to filter
    merge_distance = 4
    indexes = sorted(to_check.keys())
    for group, inds in itertools.groupby(indexes, key=lambda x: x[:2]):
        #log.debug('processing group %(group)s' % locals())
        positions = [x[2] for x in inds]
        current = [positions.pop(0)]
        while positions:
            head = current[-1]
            consider = positions[0]
            #log.debug('head is %(head)s, considering %(consider)s' % locals())
            if abs(consider - head) <= merge_distance:
                current.append(positions.pop(0))
            else:
                #log.debug('deduping %(current)s' % locals())
                for nondup in filter_dups(group, positions):
                    outp.write(nondup)
                current = [positions.pop(0)]

        for nondup in filter_dups(group, current):
            outp.write(nondup)




