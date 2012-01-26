import sys
import itertools
import random
from collections import Counter

import pysam

import amplicon
import stats


class MidAnnotator(object):
    """ Annotate BAM file with MIDs """

    # TODO: if the BAM is in the same order as the reads, can we avoid preloading read mids
    # TODO: write the barcode quality attribute

    @classmethod
    def customize_parser(cls, parser):
        parser.add_argument('--mids', type=str, help='list of mids')
        parser.add_argument('--mids-read', type=str, help='file containing the barcode contained in each read')
        parser.add_argument('--library', type=str, help='library to use in RG header')
        parser.add_argument('--platform', type=str, help='platform to use in RG header')


    def __init__(self, args, header):

        # parse the trim file of actually read mids
        read_mids = itertools.imap(
            lambda line: line.rstrip().split(),
            file(args.mids_read)
        )

        self.read_mids = dict(((y,x) for x, y in read_mids))

        # get the expected mids
        mids =  itertools.imap(
            lambda line: line.rstrip().split(),
            file(args.mids)
        )

        self.mids = dict(mids)
        self.counts = dict([(x,0) for x in self.mids.values()])

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

        read_mid = self.read_mids[read.qname]
        RG = self.match_read(read_mid)

        read.tags = read.tags + [('RG', RG), ('BC', read_mid)]
        self.counts[RG] += 1

    def report(self):
        print 'sample reads'

        for k,v in sorted(self.counts.items()):
            print k, v



class DbrAnnotator(object):
    """ Annotate BAM file with MIDs """

    # TODO: same as MidAnnotator

    @classmethod
    def customize_parser(cls, parser):
        parser.add_argument('--counters', type=str, help='Molecular counter file for each read')

    def __init__(self, args, header):

        # parse the trim file of actually read mids
        read_mids = itertools.imap(
            lambda line: line.rstrip().split(),
            file(args.counters)
        )

        self.read_mids = dict(((y,x) for x, y in read_mids))

    def __call__(self, read):
        MC = self.read_mids[read.qname]
        read.tags = read.tags + [('MC', MC)]

    def report(self):
        pass


class AmpliconAnnotator(object):

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--amps', type=str, help='amps file')
        parser.add_argument('--id-column', type=str, help='amps file', default='id')
        parser.add_argument('--amplicon-column', type=str, help='amps file', default='amplicon')
        parser.add_argument('--trim-column', type=str, help='amps file', default='trim')
        parser.add_argument('--delimiter', type=str, help='file', default='\t')
        parser.add_argument('--offset-allowed', type=int, help='file', default=10)
        parser.add_argument('--clip', action='store_true')


    def __init__(self, args, header):
        self.stats = stats.Stats('')
        self.amplicons = amplicon.load_amplicons(args.amps, self.stats, args)

        AMS = []
        for amp in self.amplicons:
            AMS.append({
                'ID': amp.external_id,
                'AC': '%s:%s-%s' % (amp.chr, amp.start, amp.end),
                'TC': '%s:%s-%s' % (amp.chr, amp.trim_start, amp.trim_end),
                'ST': str(amp.strand)
                })

        header['EA'] = AMS

    def __call__(self, read):
        for amp in self.amplicons:
            if amp.matches(read):
                # FIXME: prevent duplicate annotation
                # FIXME: amplicon mark method
                # FIXME: short circuit
                read.tags = read.tags + [('EA', amp.external_id)]
                return read


    def report(self):
        self.stats.report(sys.stdout)


def annotate(args):
    inp = pysam.Samfile(args.input)

    header = inp.header
    annotators = []
    if args.amps:
        annotators.append(AmpliconAnnotator(args, header))
    if args.mids:
        annotators.append(MidAnnotator(args, header))
    if args.dbrs:
        annotators.append(DbrAnnotator(args, header))

    assert 'SQ' in header # http://code.google.com/p/pysam/issues/detail?id=84
    oup = pysam.Samfile(args.output, 'wb', header=header)

    for read in inp:
        for annotator in annotators:
            annotator(read)
        oup.write(read)


    for a in annotators:
        a.report()


def duplicates(args):
    inp = pysam.Samfile(args.input, "rb" )
    outp = pysam.Samfile(args.output, "wb", template=inp)

    # TODO: check sorted

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
            rg, dbr = e.opt('RG'), e.opt('DB')

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




