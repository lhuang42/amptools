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

        self.amplicons = dict([(x.external_id, x) for x in self.amplicons])

    def __call__(self, samfile, outfile):
        for r in samfile:
            EA = dict(r.tags).get('ea', None)
            if EA is not None:
                clipped = self.amplicons[EA].clip(r)
                if clipped:
                    outfile.write(r)
            else:
                outfile.write(r)


def clip(args):
    """ clip primer sequences from amptools annotated BAM """
    inp = pysam.Samfile(args.input, 'rb')
    oup = pysam.Samfile(args.output, 'wb', template=inp)

    clipper = AmpliconClipper(args)
    clipper(inp, oup)
    clipper.stats.report(sys.stdout)

