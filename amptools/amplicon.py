from __future__ import print_function
import sys
import itertools
import cigar

class Amplicon(object):
        
    def __str__(self):
        return '%s:%s-%s:%s' % (self.chr, self.start, self.end, self.strand)
    
    def __init__(self, external_id=None, chr=None, start=None, end=None, 
        strand=None, trim_start=None, trim_end=None, **args):
        self.external_id = external_id
        self.chr = chr
        self.start = start 
        self.end = end
        self.strand = strand
        self.trim_start = trim_start
        self.trim_end = trim_end
    
    def reads_from(self, samfile, offset_allowed=10):
        """ Return an iterator of reads from this amplicon in the samfile """
        reads = samfile.fetch(self.chr, self.start, self.end)
        
        def filterfunc(x):
            start_correct = abs(x.pos - self.start) < offset_allowed
            end_correct = abs(x.aend - self.end) < offset_allowed
            
            if self.strand > 0: 
                orientation_correct = not x.is_reverse
                return start_correct and orientation_correct
            elif self.strand < 0: 
                orientation_correct = x.is_reverse
                return end_correct and orientation_correct
            else: 
                return start_correct or orientation_correct
                        
        return itertools.ifilter(filterfunc, reads)
    
    def pileup_dict_at_position(self, samfile, position):
        """ create a dictionary of read accession to pileup at a position """
        pileup_dict = {}
        
        for pu_column in samfile.pileup(self.chr, position, position+1):
            if pu_column.pos != position:
                continue
            
            for pu in pu_column.pileups:
                pileup_dict[pu.alignment.qname] = pu.qpos
                
        return pileup_dict
                                                                                            
    def pileup_dict_at_start(self, samfile):
        """ return a pileup at the trim start"""
        if not self.trim_start: return {}
        return self.pileup_dict_at_position(samfile, self.trim_start)
    
    def pileup_dict_at_end(self, samfile):
        if not self.trim_end: return {}
        return self.pileup_dict_at_position(samfile, self.trim_end)
    
    def clip(self, read):
        """ trim a read """
        
        first_base_pos = self.start_trim_dict.get(read.qname, 0)
        last_base_pos = self.end_trim_dict.get(read.qname, -1)
    
        # cache original details, otherwise they go to None
        seq, qual, cig, end_pos = read.seq, read.qual, read.cigar, read.aend
        
        if first_base_pos: 
            read.seq = seq[first_base_pos:]
            read.qual = qual[first_base_pos:]            
            read.cigar = cigar.trim_cigar(cig, len(read.seq), start=True)
            read.pos = end_pos - cigar.ref_length(read.cigar)
            assert read.aend == end_pos
            
            # recache new stuff and recalculate last base
            seq, qual, cig, end_pos = read.seq, read.qual, read.cigar, read.aend
            
            if last_base_pos != -1:
                last_base_pos = last_base_pos - first_base_pos
            
        if last_base_pos != -1: 
            read.seq = seq[:last_base_pos]
            read.qual = qual[:last_base_pos]
            read.cigar = cigar.trim_cigar(cig, len(read.seq))

        return read
    
    
    def clipped_reads(self, samfile):
        """ return all the reads of this amplicon clipped """

        reads = self.reads_from(samfile)
        self.start_trim_dict = self.pileup_dict_at_start(samfile)
        self.end_trim_dict = self.pileup_dict_at_end(samfile)
        
        return itertools.imap(self.clip, reads)
