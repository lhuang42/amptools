"""
Methods for handling cigar lists.

Pysam returns a list of [operation, bases] representing the alignment of a
read to a reference
"""
from __future__ import print_function
import sys
import logging; log = logging.getLogger(__name__)

# cigar operations
MATCH = 0 # can be match or substitution
INS = 1 # insert to reference
DEL = 2 # deletion to reference
SKIP = 3
SOFT_CLIP = 4
HARD_CLIP = 5
PADDING = 6

# indexes of operation or number of bases
OP = 0
BASES = 1


def read_length(cigar):
    """ Return the number of read bases represented by cigar """
    return sum(
        [
            x[BASES] for x in cigar
            if x[OP] in [MATCH, INS, SOFT_CLIP]

        ]
    )

def ref_length(cigar):
    """ Return the number of reference bases represented by the cigar """
    return sum(
        [x[BASES] for x in cigar
        if x[OP] in [MATCH, DEL]]
    )

def remove_soft(cigar, seq, qual):
    """ Remove soft clipped bases from a cigar string and sequence """
    if cigar[0][OP] == SOFT_CLIP:
        sc = cigar.pop(0)
        seq = seq[sc[BASES]:]
        qual = qual[sc[BASES]:]
    if cigar[-1][OP] == SOFT_CLIP:
        sc = cigar.pop()
        seq = seq[:-sc[BASES]]
        qual = qual[:-sc[BASES]]

    return cigar, seq, qual

def trim_cigar(cigar, n, start=False):
    """ Trim cigar until it represents n read bases, inserting hard clips

        Defaults to trimming at the 3' end of the read, set start=True to
        trim from the 5' end (start) of the read
    """
    # if we need to trim from the start, reverse the cigar, trim and reverse the result
    if start:
        cigar = list(reversed(cigar))
        cigar = trim_cigar(cigar, n)
        return list(reversed(cigar))

    # check we can handle these ops
    all_ops = [x[OP] for x in cigar]
    if PADDING in all_ops or SKIP in all_ops:
        raise NotImplementedError

    to_trim = read_length(cigar) - n
    if to_trim == 0:
        return cigar

    # cache the clips
    start_clip, end_clip = [], [(HARD_CLIP, 0)]

    if cigar[0][OP] == HARD_CLIP:
        start_clip = [cigar.pop(0)]
    if cigar[-1][OP] ==  HARD_CLIP:
        end_clip = [cigar.pop()]

    assert to_trim > 0
    end_clip = [(HARD_CLIP, end_clip[0][BASES] + to_trim)]

    while to_trim:

        op = cigar[-1][OP]
        bases = cigar[-1][BASES]

        if op in [MATCH, INS, SOFT_CLIP]:

            # not enough bases
            if bases - to_trim < 0:
                to_trim = to_trim - bases
                cigar.pop()

            # exact number of bases
            elif bases - to_trim == 0:
                to_trim = 0
                cigar.pop()

            # too many bases
            else:
                # print 'trimmed', to_trim
                cigar[-1] = (op, bases - to_trim)
                to_trim = 0


        elif op == DEL:
            cigar.pop()

        else:
            raise Exception('bad cigar element in trimming: %s' % cigar[-1])

    # remove final deletions as they will affect the placement of
    # reversed reads
    while cigar[-1][OP] == DEL:
        cigar.pop()

    assert read_length(cigar) == n, '%s is not length %s' % (cigar, n)
    return start_clip + cigar + end_clip
