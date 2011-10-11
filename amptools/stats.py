from __future__ import print_function, division

from collections import Counter
import csv

class Stats(object):
    
    def __init__(self, cmdline): 
        self.cmdline = cmdline
        self._start_trims = Counter()
        self._end_trims = Counter()
        self._matches = Counter()
        self.reads = 0
        
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
        for eid in sorted(set(self._start_trims.keys() + self._end_trims.keys() + self._matches.keys() )): 
            writer.writerow([eid, 
                str(self._matches[eid]), 
                str(self._start_trims[eid]), 
                str(self._end_trims[eid])])
            
        reads = self.reads
        matched = sum(self._matches.values())
        matched_pc = 100 * matched / reads
        print(file=stream)
        print("matched %(matched)s/%(reads)s %(matched_pc).2f%% of alignments" % locals(), file=stream)
    