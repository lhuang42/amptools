from __future__ import print_function, division

from collections import Counter
import csv
from rpy2 import robjects

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
    gtp = 0.05 + (0.9*gt/2)
    
    
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
    
    ra = robjects.IntVector([x['RA'] for x in calls])
    aa = robjects.IntVector([x['AA'] for x in calls])
    gt = robjects.IntVector(map(sum, [x['GT'] for x in calls]))
    test_val, ab = ll_test(ra, aa, gt)
    # print('test value -ve prefers gt, +ve prefers error', test_val, test_val < 0)
    return test_val < 0, test_val, ab
    