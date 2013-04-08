import sys
import itertools
import logging; log = logging.getLogger(__name__)

import json
from collections import Counter

import pysam
import ngram

import amplicon
import stats
import subprocess

from Bio import AlignIO

TAG_COUNT = 'mc'
TAG_AMP = 'ea'
TAG_BC = 'BC'
TAG_RG = 'RG'


def get_alignment(reads):
    ind = range(0, len(reads))
    i = zip(ind, reads)
    fa = ['>' + str(x) + '\n' + y.seq + '\n' for x,y in i]

    p = subprocess.Popen(['muscle', '-in', '-', '-out', '-', '-quiet',
                          '-clwstrict'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    for seq in fa:
        p.stdin.write(str(seq))
    p.stdin.close()

    return AlignIO.read(p.stdout, 'clustal')

def get_mode_bq(reads_in):
    try:
        qlen = max([len(x.seq) for x in reads_in])
    except ValueError:
        return ''

    bqs = [zip(range(0, len(x.qual)), list(x.qual)) for x in reads_in]
    
    quals = {}
    for x in range(0, qlen):
        quals[x] = ''

    for read in bqs:
        for pos, bq in read:
            quals[pos] += bq

    for k,v in quals.items():
        quals[k] = Counter(v).most_common(1)[0]
        if len(quals[k]) == 0:
            quals[k] = ''
        else:
            quals[k] = quals[k][0]
            
    return ''.join([quals[x] for x in quals])
    

def consensus(reads, THRESHOLD=0.80):
    r1con, r2con = None, None
    con_r1, con_r2 = None, None

    r1 = [rx for rx in reads if rx.is_read1]
    r2 = [rx for rx in reads if rx.is_read2]

    rule = lambda xx: xx.most_common(1)[0][1] / float(sum(xx.values())) > THRESHOLD

    if len(r1) > 1:
        r1align = get_alignment(r1)
        r1composition = [Counter(r1align[:,pos]) for pos in range(r1align.get_alignment_length())]
        
        for x in r1composition:
            del x['-']
        r1con = ''.join([x.most_common(1)[0][0] if rule(x) else 'X' for x in r1composition])
        
    if len(r2) > 1:
        r2align = get_alignment(r2)
        r2composition = [Counter(r2align[:,pos]) for pos in range(r2align.get_alignment_length())]
        
        for x in r2composition:
            del x['-']
        r2con = ''.join([x.most_common(1)[0][0] if rule(x) else 'X' for x in r2composition])

    r1bq = get_mode_bq(r1)
    r2bq = get_mode_bq(r2)

    if r1con:
        if len(r1bq) < len(r1con):
            xbq = len(r1con) - len(r1bq)
            r1bq = r1bq + ('B' * xbq)
            
        con_r1 = pysam.AlignedRead()
        con_r1.is_duplicate = False
        con_r1.seq = r1con
        con_r1.flag = 99
        con_r1.rname = r1[0].rname
        con_r1.pos = min([x.pos for x in r1])
        con_r1.mapq = max([x.mapq for x in r1])
        con_r1.mrnm = r1[0].mrnm
        con_r1.qual = r1bq
        con_r1.tags = r1[0].tags + [('NR',len(r1))]
        con_r1.cigar = [(0,len(con_r1.seq))]
        con_r1.qname = 'con_' + str(con_r1.pos + 1) + '_' + con_r1.opt('mc') + '_' + con_r1.opt('RG') + ':R1'
        con_r1.is_read1 = True
        con_r1.mate_is_reverse = True
        con_r1.is_proper_pair = True
        del r1[-1]
        for x in r1:
            x.is_duplicate = True
    else:
        try:
            r1[0].is_duplicate = False
        except IndexError:
            log.debug ('No read 1')

    if r2con:
        if len(r2bq) < len(r2con):
            xbq = len(r2con) - len(r2bq)
            r2bq = ('B' * xbq) + r2bq
        con_r2 = pysam.AlignedRead()
        con_r2.is_duplicate = False
        con_r2.seq = r2con
        con_r2.flag = 147
        con_r2.rname = r2[0].rname
        con_r2.pos = min([x.pos for x in r2]) - 1
        con_r2.mapq = max([x.mapq for x in r2])
        con_r2.mrnm = max([x.mrnm for x in r2])
        con_r2.qual = r2bq
        con_r2.tags = r2[0].tags + [('NR',len(r2))]
        con_r2.cigar = [(0,len(con_r2.seq))]
        con_r2.qname = 'con_' + str(con_r2.pos) + '_' + con_r2.opt('mc') + '_' + con_r2.opt('RG') + ':R2'
        con_r2.is_reverse = True
        con_r2.is_read2 = True
        con_r2.mate_is_reverse = False
        con_r2.is_proper_pair = True
        del r2[-1]
        for x in r2:
            x.is_duplicate = True
    else:
        try:
            r2[0].is_duplicate = False
        except IndexError:
            log.debug('No read 2')

    # fix mate
    if con_r1 and con_r2:
        con_r1.mpos = con_r2.pos
        con_r2.mpos = con_r1.pos
        con_r1.isize = (con_r2.pos + len(con_r2.seq)) - con_r1.pos
        con_r2.isize = con_r1.isize

    reads = [] + r1 + r2
    if con_r1:
        reads.append(con_r1)
    if con_r2:
        reads.append(con_r2)

    return reads
    

    

def duplicates(args):
    """ Mark duplicates using a molecular counter.

        The file should contain MC tags.  Duplicates are detected by looking
        for reads in the same start position and orientation and with the same MC tag.

        WARNING: this unsorts your input
    """
    inp = pysam.Samfile(args.input, "rb" )
    outp = pysam.Samfile(args.output, "wb", template=inp)

    # TODO: check sorted
    # TODO: report duplicates/unannotated stats
    # TODO: add sort option when output set

    to_check = {}

    def filter_dups(mid, counter, amplicon):
        """ filter a set of start positions for the best read by mapping quality """
        entries = to_check.pop((mid, counter, amplicon))
        nr = len(entries)
        log.debug('using %d entries from %s, %s, %s' % (nr, mid, counter, amplicon))

        # return the best read in each group
        for dup in consensus(entries):
            yield dup


    for (i, entry) in enumerate(inp):
        if (i % 100) == 0:
            log.info('processed %(i)s reads' % locals())

        # append group to the to_check index
        position_index = (entry.opt('RG'), entry.opt('mc'), entry.opt('ea'))
        to_check[position_index] = to_check.get(position_index, []) + [entry]

    indexes = sorted(to_check.keys())
    log.debug(indexes)        
    for mid, counter, amplicon in indexes:
        log.debug('deduping %(mid)s, %(counter)s, %(amplicon)s' % locals())
        for nondup in filter_dups(mid, counter, amplicon):
            outp.write(nondup)

