from __future__ import print_function, division
import logging; log = logging.getLogger(__name__)

from collections import Counter
import csv
import sys
from rpy2 import robjects
from rpy2.robjects.packages import importr
import amplicon

import pysam

class Stats(object):

    def __init__(self, cmdline):
        self.cmdline = cmdline
        self._start_trims = Counter()
        self._end_trims = Counter()
        self._matches = Counter()
        self.reads = 0
        self.eids = []

    def start_trim(self, eid):
        self._start_trims[eid] += 1

    def end_trim(self, eid):
        self._end_trims[eid] += 1

    def match(self, eid):
        self._matches[eid] += 1

    def report(self, stream):

        print('amptools version xx', file=stream)
        print('invocation:', self.cmdline, file=stream)
        print(file=stream)
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow('amplicon matches start_trims end_trims'.split())
        for eid in sorted(self.eids):
            writer.writerow([eid,
                str(self._matches[eid]),
                str(self._start_trims[eid]),
                str(self._end_trims[eid])])

        reads = self.reads
        if reads:
            matched = sum(self._matches.values())
            matched_pc = 100 * matched / reads
            print(file=stream)
            print("matched %(matched)s/%(reads)s %(matched_pc).2f%% of alignments" % locals(), file=stream)



ll_test = robjects.r('''
function(ra, aa, gt, diag=F) {
    ra_sum = sum(ra)
    aa_sum = sum(aa)
    ab = aa_sum / (ra_sum + aa_sum)
    gtp = 0.5 + (0.48*(gt-1))

    error_likelihood = log(dbinom(aa, ra+aa, ab))
    gt_likelihood = log(dbinom(aa, ra+aa, gtp))

    if (diag) {
        print(ra)
        print(aa)
        print(gtp)
        print(ab)
        print(error_likelihood)
        print(gt_likelihood)
    }
    error_likelihood = sum(error_likelihood)
    gt_likelihood = sum(gt_likelihood)
    if (F) {
        print(error_likelihood)
        print(gt_likelihood)

        print(exp(error_likelihood))
        print(exp(gt_likelihood))
    }
    c(error_likelihood - gt_likelihood, ab)
}
''')

def bias_test(calls):


    calls = [x for x in calls if x.called]

    #TODO: single genotype assumption
    try:
        # freebayes
        ra = robjects.IntVector([x['RO'][0] for x in calls])
        aa = robjects.IntVector([x['AO'][0] for x in calls])
    except KeyError:
        # GATK
        ra = robjects.IntVector([x['AD'][0] for x in calls ])
        aa = robjects.IntVector([x['AD'][1] for x in calls ])


    gt = robjects.IntVector([x.gt_type for x in calls])
    # print(ra, aa, gt)
    test_val, ab = ll_test(ra, aa, gt)
    #print('test value -ve prefers gt, +ve prefers error', test_val, test_val < 0, ab)

    return test_val < 0, test_val, ab

def neg_binom_fit(reads):
    robjects.r.options(warn=-1)
    reads = robjects.IntVector(reads)
    mass = importr('MASS')
    fit = mass.fitdistr(reads, 'negative binomial')
    print('Negative binomial fit:')
    print(fit)
    qvals = robjects.FloatVector([0,0.01,0.05,0.1,0.25,0.5])
    print('Quantiles:')
    print(robjects.r.quantile(reads, qvals))


def _cluster_center(calls):
    if not calls: return None
    ra = sum([x['RO'] for x in calls])
    aa = sum([x['AO'] for x in calls])
    return ra / (ra + aa)

def _mean_distance_to_center(calls):
    if not calls: return None
    center = _cluster_center(calls)
    for x in calls:
        if x['RO'] + x['AO'] == 0:
            print ('x'*80, '\n', x)
    posns = [(x['RO']/( x['RO'] + x['AO'] )) for x in calls if (x['RO'] + x['AO']) > 0 ]
    dists = [abs(x - center) for x in posns]
    return sum(dists) / len(dists)

def cluster_separation(calls):

    hom = [x for x in calls if x['GT'] == (0,0)]
    het = [x for x in calls if x['GT'] == (0,1)]
    var = [x for x in calls if x['GT'] == (1,1)]

    # print('dist', len(hom), len(het), len(var))

    values = [0]
    for (left, right) in [(hom, het), (het, var)]:
        left_center, right_center = _cluster_center(left), _cluster_center(right)
        if left_center is None or right_center is None:
            continue
        dist_between_means = abs(left_center - right_center)
        score = _mean_distance_to_center(left) + _mean_distance_to_center(right)
        score = score / dist_between_means
        values.append(score)
    # print('vals', values)

    return max(values)


amp_bias_test = ('''
function(AOs, DPs) {which.max(c(dbinom(AOs, DPs, 0.01),  dbinom(AOs, DPs, 0.5), dbinom(AOs, DPs, 0.99))) - 1}
''')

def check_amplicon_bias(calls, counts, amps, amp_counter):
    """Iterate through a number of VCF calls and write AMPC into each record in place.
        AMPC = supporting amplion count, i.e. the number of amplicons where
        the call for the amplicon matches the overall call.
     """

    # this R callout needs to be vectorized
    for call in calls:
        if call['_GT'] is None:
            continue

        gt = sum(call['_GT'])

        supporting = 0
        ampcount = 0
        for amp in amps:
            ro = counts.get((call['name'], amp, True), 0)
            ao = counts.get((call['name'], amp, False), 0)
            # print ('test', call['name'], amp, ao, ao+ro, amp_bias_test(ao, ao+ro)[0])

            if ro + ao > 0:
                ampcount += 1
                if amp_bias_test(ao, ao+ro)[0] == gt:
                    supporting += 1
                else:
                    amp_counter[amp] += 1


        call['AMPS'] = supporting
        call['AMPC'] = ampcount




def coverage(args):
    counts = {}
    inp = pysam.Samfile(args.input)

    uniq = {}
    reads = {}

    total = 0
    stats = Stats('')
    amplicons = amplicon.load_amplicons_from_header(inp.header, stats, None)
    libs = {}

    for rg in inp.header['RG']:
        for amp in amplicons:
            key = rg['ID'], amp.external_id
            reads[key] = uniq[key] = 0
            libs[rg['ID']] = rg.get('LB', None)

    for r in inp:
        total += 1
        tags = dict(r.tags)
        try:
            key = tags['RG'], tags['ea']
        except KeyError:
            continue
        try:
            reads[key] += 1
            if not r.is_duplicate:
                uniq[key] += 1
        except KeyError:
            logging.debug('unexpected key')

    total_ot = sum(reads.values())
    total_uniq = sum(uniq.values())

    total_ot_p = (100.0 * total_ot) / total

    try:
        reads_per_counter = float(total_ot) / total_uniq
    except ZeroDivisionError:
        reads_per_counter = 0

    print('total %(total)s reads, on target %(total_ot)s, uniq %(total_uniq)s' % locals(), file=sys.stderr)
    print('on target %3.2f%%' % total_ot_p, file=sys.stderr)
    print('on target reads per counter: %2.2f' % reads_per_counter, file=sys.stderr)

    if args.control:
        total_control = sum([reads[x] for x in reads if x[0] == args.control])
        control_p = (100*total_control)/total
        print('control reads %(total_control)s, %(control_p)f%%' % locals(), file=sys.stderr)

    if not args.output:
        out = csv.writer(sys.stdout)
    else:
        out = csv.writer(open(args.output, 'w'))
    out.writerow(['rg', 'lib', 'amp', 'unique', 'reads'])
    for (rg, amp) in sorted(reads):
        if rg != args.control:
            key = (rg, amp)
            out.writerow(map(str, (rg, libs[rg], amp, uniq[key], reads[key])))

