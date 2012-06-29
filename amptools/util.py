from collections import namedtuple
import vcf.filters
import pysam


PILEUP_DEPTH = 300000


class AmpliconFilter(vcf.filters.Base):

    name = 'ac'
    description = "Test for at least two amplicons for a novel SNP"

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--reads', type=str,
                help='Samfile')

    def __init__(self, args):
        self.reads = pysam.Samfile(args.reads)
        self.threshold = 10

    def __call__(self, entry):
        if entry.is_monomorphic or entry.ID:
            return None
        if entry.is_indel:
            raise NotImplementedError()

        counts = {}
        ab = {}
        posn = entry.POS - 1
        for puc in self.reads.pileup(entry.CHROM, posn-1, posn+1, max_depth=PILEUP_DEPTH):
            if puc.pos != posn:
                continue

            for pu in puc.pileups:
                if pu.is_del:
                    continue
                base = pu.alignment.query[pu.qpos]
                ea = dict(pu.alignment.tags).get('ea', None)
                ab[base] = ab.get(base, 0) + 1
                if ea and base == entry.ALT[0]:
                    counts[ea] = counts.get(ea, 0) + 1
        amps = 0
        for amp, count in counts.items():
            if count > self.threshold:
                amps += 1
        if amps < 2:
            return True
