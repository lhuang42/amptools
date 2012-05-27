from collections import Counter, namedtuple

import vcf
import vcf.filters

import pysam

import stats

Obs = namedtuple('Observation', ['rg', 'amp', 'is_ref'])

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




def _dep(entry, samfile):
    # TODO: single alternate allele assumption
    if len(entry.REF) != len(entry.ALT[0]):
        is_indel = True

        if len(entry.REF) > len(entry.ALT[0]):
            is_del = True
        else:
            is_del = False

    else:
        is_del = False
        is_indel = False
        mnp_length = len(entry.REF)

    counts = Counter()
    read_meta = {}
    read_bases = {}

    posns = range(entry.POS - 1, entry.POS - 1  + len(entry.REF))



    for puc in samfile.pileup(entry.CHROM, posns[0], posns[-1]+1, max_depth=PILEUP_DEPTH):
        if puc.pos in posns:

            for pu in puc.pileups:
                acc = pu.alignment.qname
                # if acc not in ['HAV2UQ001B2AB7', 'HAV2UQ001BYY8P']: continue

                if pu.alignment.qname not in read_meta:
                    tags = dict(pu.alignment.tags)
                    read_meta[acc] = tags.get('RG', None), tags.get('AM', None)

                if pu.indel <= 0 :
                    base = pu.alignment.query[pu.qpos]
                elif pu.indel > 0 and not pu.is_del:
                    # print 'woot'
                    base = pu.alignment.query[pu.qpos:pu.qpos+pu.indel+1]

                if pu.is_del:
                    base = ''
                read_bases[acc] = read_bases.get(acc, '') + base

                # if pu.indel or pu.is_del:
                #     print acc,puc.pos, base, pu.qpos, pu.indel, pu.is_del, pu.alignment.query[pu.qpos]


    # print read_bases.items()[:5]

    for acc in read_bases:
        read = read_bases[acc]
        rg, amp = read_meta[acc]

        # if read != entry.REF:
        #     print acc, read, read in (entry.REF, entry.ALT[0])

        if read in (entry.REF, entry.ALT[0]):
            observation = Obs(rg=rg, amp=amp, is_ref=read==entry.REF)
            counts[observation] += 1

    return counts



class SeqErrFilter(vcf.filters.Base):

    name = 'eb'
    description = "Test for error bias"

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--errlr', type=int, default=-10,
                help='Filter sites above this error log odds ratio')

    def __init__(self, args):
        self.threshold = args.errlr

    def __call__(self, record):
        if record.is_monomorphic:
            return None
        passed, tv, ab = stats.bias_test(record.samples)
        if tv > self.threshold:
            return tv



