import sys
import itertools
from collections import Counter

import pysam

import amplicon
import stats


class MidAnnotator(object):
    """ Annotate BAM file with MIDs """

    @classmethod
    def customize_parser(cls, parser):
        parser.add_argument('--mids', type=str, help='mids file')
        parser.add_argument('--trim', type=str, help='trim file')


    def __init__(self, args):

        # parse the trim file of actually read mids
        read_mids = itertools.imap(
            lambda line: line.rstrip().split(),
            file(args.trim)
        )

        self.read_mids = dict(((y,x) for x, y in read_mids))

        # get the expected mids
        mids =  itertools.imap(
            lambda line: line.rstrip().split(),
            file(args.mids)
        )

        self.mids = dict(mids)
        self.counts = dict([(x,0) for x in self.mids.values()])

    def match_read(self, read_mid):
        return self.mids[read_mid]

    def __call__(self, read):

        read_mid = self.read_mids[read.qname]
        RG = self.match_read(read_mid)

        read.tags = read.tags + [('RG', RG), ('MI', read_mid)]
        self.counts[RG] += 1

    def report(self):
        print 'sample reads'

        for k,v in sorted(self.counts.items()):
            print k, v



class DbrAnnotator(object):
    """ Annotate BAM file with MIDs """

    @classmethod
    def customize_parser(cls, parser):
        parser.add_argument('--dbrs', type=str, help='dbr file')

    def __init__(self, args):

        # parse the trim file of actually read mids
        read_mids = itertools.imap(
            lambda line: line.rstrip().split(),
            file(args.dbrs)
        )

        self.read_mids = dict(((y,x) for x, y in read_mids))

    def __call__(self, read):
        DB = self.read_mids[read.qname]
        read.tags = read.tags + [('DB', DB)]

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


    def __init__(self, args):
        self.stats = stats.Stats('')
        self.amplicons = amplicon.load_amplicons(args.amps, self.stats, args)

    def __call__(self, read):
        for amp in self.amplicons:
            if amp.matches(read):
                read.tags = read.tags + [('AM', amp.external_id)]
                return read


    def report(self):
        self.stats.report(sys.stdout)


def annotate(args):
    inp = pysam.Samfile(args.input, 'rb')
    oup = pysam.Samfile(args.output, 'wb', template=inp)

    annotators = []
    if args.amps:
        annotators.append(AmpliconAnnotator(args))
    if args.mids:
        annotators.append(MidAnnotator(args))
    if args.dbrs:
        annotators.append(DbrAnnotator(args))

    for read in inp:
        for annotator in annotators:
            annotator(read)
        oup.write(read)


    for a in annotators:
        a.report()


