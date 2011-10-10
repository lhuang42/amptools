from __future__ import print_function

import sys
import cmdln
import pysam
from amplicon import Amplicon
from ucsc import Interval    

class Amptools(cmdln.Cmdln):
    name = "amptools"


    @cmdln.option("-o", "--outfile", action="store", default='-',
                  help="Output file (default stdout)")
    
    def do_clip_one(self, subcmd, opts, bamfile, amplicon_location, trim_location):
        samfile = pysam.Samfile(bamfile, 'rb')
        outfile = pysam.Samfile(opts.outfile, 'wb', template=samfile)

        amp_loc = Interval.from_string(amplicon_location)
        trim_loc = Interval.from_string(trim_location)
        
        if not amp_loc.contains(trim_loc):
            print('trim location not contained in amplicon location, impossible trim', file=sys.stderr)
            sys.exit(1)
        
        amplicon = Amplicon(chr=amp_loc.chrom, start=amp_loc.start, end=amp_loc.end,
            trim_start = trim_loc.start, trim_end=trim_loc.end
        )
        
        trimmed = amplicon.clipped_reads(samfile)
        map(outfile.write, trimmed)