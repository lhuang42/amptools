import sys

import pysam

import stats
import amplicon

class AmpliconClipper(object):

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
        self.samfile = pysam.Samfile(args.input)
        self.amplicons = amplicon.load_amplicons_from_header(self.samfile.header, self.stats, self.samfile)

    def __call__(self, samfile, outfile):
        clipped = {}
        for amplicon in self.amplicons:
            trimmed = amplicon.clipped_reads(samfile)
            for t in trimmed:
                outfile.write(t)
                clipped[t.qname] = True

        for r in samfile.fetch():
            if r.qname not in clipped:
                outfile.write(r)


def clip(args):
    """ clip primer sequences from amptools annotated BAM """
    inp = pysam.Samfile(args.input, 'rb')
    oup = pysam.Samfile(args.output, 'wb', template=inp)

    clipper = AmpliconClipper(args)
    clipper(inp, oup)
    clipper.stats.report(sys.stdout)

