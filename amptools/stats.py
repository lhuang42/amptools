from __future__ import print_function

from collections import defaultdict
import csv

class Stats(object):
    
    def __init__(self, cmdline): 
        self.cmdline = cmdline
        self._start_trims = defaultdict(int)
        self._end_trims = defaultdict(int)
        
    def start_trim(self, eid):
        self._start_trims[eid] += 1
        
    def end_trim(self, eid):
        self._end_trims[eid] += 1
        
    def report(self, stream):
        
        print('amptools version xx', file=stream)
        print('invocation:', self.cmdline, file=stream)
        print(file=stream)
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow('amplicon start_trims end_trims'.split())
        for eid in sorted(set(self._start_trims.keys() + self._end_trims.keys())): 
            writer.writerow([eid, str(self._start_trims[eid]), str(self._end_trims[eid])])
            
    