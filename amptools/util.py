from collections import Counter, namedtuple

import vcf

import stats

Obs = namedtuple('Observation', ['rg', 'amp', 'is_ref'])

PILEUP_DEPTH = 300000

def amplicon_distribution(entry, samfile):
    """ Compute the distribution of a VCF entry over different amplicons

        Given a vcf entry and a samfile, this method looks at each read and
        decides whether it consitutes a reference of alternate allele.
        It returns a Counter containing (rg, amplicon, is_ref)
    """

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



class SeqErrFilter(vcf.Filter):

    name = 'seq_err_filter'
    short_name = 'selr'
    description = "Test probability site is a sequencing error using constant rate vs genotypes"

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--errlr', type=int, default=-10,
                help='Filter sites above this error log odds ratio')

    def __init__(self, args):
        self.threshold = args.errlr

    def __call__(self, record):
        passed, tv, ab = stats.bias_test(record.samples)
        if tv > self.threshold:
            return tv



def apply_filters(record, error_bias=-10, genotype_quality = 30, min_ampc=None):

    # cannot handle composite classes
    if len(record.ALT) > 1:
        record.add_filter('composite')
        return

    # biased site filter: tests for all sites coming from an error probability
    if error_bias:
        passed, tv, ab = stats.bias_test(record.samples)
        record.add_info('ERRLR', tv)#'=%(tv).2f;ERRAB=%(ab).2f' % locals())
        if tv >= error_bias:
            record.add_filter('ERRLR=%(tv).2f' % locals())
            # return

    if min_ampc:
        # check for a variant that is only supported by a single amplicon
        # we do not allow novel variants only seen by one amplicon
        ampcs = [int(x['AMPC']) for x in record.samples if x['_GT'] is not None and sum(x['_GT']) != 0]
        if all([x<2 for x in ampcs]):
            if record.ID == '.':
                record.add_filter('AMPC=%s' % min(ampcs))

        # check for amplicon support for each call
        # for problem calls, set the GQ to 0
        for call in record.samples:
            if call['_GT'] is None or call['_GT'] == 0:
                continue

            if call['AMPS'] <  call['AMPC']:
                call['GQ'] = 0

        # if the
        # ampss = [int(x['AMPS']) for x in record.samples if x['_GT'] is not None and sum(x['_GT']) != 0]
        # if not any([x>=min_ampc for x in ampss]):
        #     record.add_filter('AMPS=%s' % max(ampss))



    if genotype_quality:
        quals = [x['GQ'] for x in record.samples if x['_GT'] is not None and sum(x['_GT']) != 0]
        if not any([x>=genotype_quality for x in quals]):
            record.add_filter('max_GQ=%s' % max(quals))
            # return


    #
    # if homopolymer_indel:
    #     is_ins, is_del = is_homo_ins(fields), is_homo_del(fields)
    #     line = add_info(line, 'HOMOINDEL=%s' % ("T" if is_ins or is_del else 'F'))
    #
    #     if 'REPEAT' in fields.INFO:
    #         repeats = parse_repeats(fields.INFO['REPEAT'])
    #
    #         if is_del:
    #             ref = fields.REF
    #             if ref == ref[0] * len(ref):
    #                 ref = ref[0]
    #             ref_count = repeats.get(ref, 0)
    #         else:
    #             ref_count = repeats.get(fields.REF, 0)
    #         # print is_ins, is_del, ref_count
    #         if (is_ins or is_del) and ref_count > max_indel_repeat:
    #             line = add_info(line, 'REPEAT=%s' % ref_count)



