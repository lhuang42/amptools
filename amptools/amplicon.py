from __future__ import print_function
import sys
import itertools
import csv
import cigar
import json
from fastinterval import Interval
import logging; log = logging.getLogger(__name__)
# set the pileup engine to allow 1500 samples at depth of 200
PILEUP_MAX_DEPTH = 200 * 1500



def load_amplicons(design, stats, opts):
    amplicons = []
    for row in csv.DictReader(file(design, 'U'), delimiter=opts.delimiter):
        amp_loc = Interval.from_string(row[opts.amplicon_column])
        trim_loc = Interval.from_string(row[opts.trim_column])


        if not trim_loc in amp_loc:
            print('trim location not contained in amplicon location, impossible trim', file=sys.stderr)
            sys.exit(1)

        amplicon = Amplicon(
            chr=amp_loc.chrom, start=amp_loc.start, end=amp_loc.end, strand=amp_loc.strand,
            trim_start = trim_loc.start, trim_end=trim_loc.end,
            external_id = row[opts.id_column], stats = stats,
            offset_allowed=opts.offset_allowed,
        )
        amplicons.append(amplicon)

    # TODO: check that amplicons can be uniquely resolved

    return amplicons


def load_amplicons_from_header(header, stats, samfile, clip=True, load_pileups=True):
    amplicons = []
    for row in header['CO']:
        try:
            row = json.loads(row)
        except:
            continue
        if row.get('type', None) != 'ea':
            continue

        amp_loc = Interval.from_string(row['ac'])
        trim_loc = Interval.from_string(row['tc'])
        strand = row.get('st')
        strand = int(strand) if strand != 'None' else 0

        amplicon = Amplicon(
            chr=amp_loc.chrom, start=amp_loc.start, end=amp_loc.end, strand=strand,
            trim_start = trim_loc.start, trim_end=trim_loc.end,
            external_id = row['id'], stats = stats,
            offset_allowed=10,
        )
        amplicons.append(amplicon)

    return amplicons


class Amplicon(object):

    def __str__(self):
        return '%s:%s-%s:%s' % (self.chr, self.start, self.end, self.strand)

    def __init__(self, external_id=None, chr=None, start=None, end=None,
        strand=None, trim_start=None, trim_end=None, stats=None,
        offset_allowed=10):
        self.external_id = external_id
        self.chr = chr
        self.start = start
        self.end = end
        self.strand = strand
        self.trim_start = trim_start
        self.trim_end = trim_end
        self.stats = stats
        self.offset_allowed = offset_allowed
        self.stats.eids.append(self.external_id)

    # TODO: use EA tag here?
    def reads_from(self, samfile):
        """ Return an iterator of reads from this amplicon in the samfile """
        reads = samfile.fetch(self.chr, self.start, self.end)
        return itertools.ifilter(self.matches, reads)

    def matches(self, read):
        """ return True if the read looks like it came from this amplicon"""
        if read.is_unmapped:
            return False

        start_correct = abs(read.pos - self.start) < self.offset_allowed
        end_correct = abs(read.aend - self.end) < self.offset_allowed

        if self.strand and self.strand > 0:
            orientation_correct = not read.is_reverse
            match =  start_correct and orientation_correct
        elif self.strand and self.strand < 0:
            orientation_correct = read.is_reverse
            match = end_correct and orientation_correct
        else:
            if read.is_reverse:
                match = end_correct
            else:
                match = start_correct

        if match:
            self.stats.match(self.external_id)

        return match


    def _find_position(self, pos, posns, lower=True):

        # find the pairing to a given position in the reference
        for x in posns:
            if x[1] == pos:
                match = x
                break
        else:
            return None

        # if this is aligned to a read base, return that
        if match[0] is not None:
            return match[0]

        # otherwise, walk the alignment to find the next base
        index = posns.index(match)
        while match[0] is None:

            if lower:
                index += 1
                if index >= len(posns):
                    return None
            else:
                index -=1
                if index < 0:
                    return None

            match = posns[index]
        return match[0]


    def clip(self, read):
        """ clip a read to the trim locations of this amplicon """

        # calculate where to clip
        posns = read.positions
        first_base_pos = self._find_position(self.trim_start, posns)
        last_base_pos = self._find_position(self.trim_end, posns, lower=False)

        # exlude primer only reads
        # TODO: may be better to check the read start end for this
        if not (first_base_pos or last_base_pos):
            return None

        # cache original details, otherwise they go to None
        seq, qual, cig, end_pos = read.seq, read.qual, read.cigar, read.aend

        # pileup position ignores soft clips, so remove
        cig, seq, qual = cigar.remove_soft(cig, seq, qual)

        if first_base_pos:
            self.stats.start_trim(self.external_id)
            read.seq = seq[first_base_pos:]
            if qual:
                read.qual = qual[first_base_pos:]
            read.cigar = cigar.trim_cigar(cig, len(read.seq), start=True)
            read.pos = end_pos - cigar.ref_length(read.cigar)
            assert read.aend == end_pos

            # recache new stuff and recalculate last base
            seq, qual, cig, end_pos = read.seq, read.qual, read.cigar, read.aend

            if last_base_pos is not None:
                last_base_pos = last_base_pos - first_base_pos

        if last_base_pos is not None:
            self.stats.end_trim(self.external_id)
            read.seq = seq[:last_base_pos]
            if qual:
                read.qual = qual[:last_base_pos]
            read.cigar = cigar.trim_cigar(cig, len(read.seq))

        return read

    def mark(self, read):
        read.tags = (read.tags or []) + [('EA', self.external_id)]
        return read

