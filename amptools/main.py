from __future__ import print_function

import sys
import csv
import cmdln
import pysam
from amplicon import Amplicon
from ucsc import Interval    
from stats import Stats

class Amptools(cmdln.Cmdln):
    """Usage:
            amptools SUBCOMMAND [ARGS...]
            amptools help SUBCOMMAND

        Amptools provides a set of tools for primer clipping SAM files to remove
        primer sites and non genomic regions in the files. Amptools requires that 
        the source file is sorted and indexed because it uses the pileup engine 
        to find trimpoints. 

        ${command_list}
        ${help_list}

        For additional information, bug reports and development, see 
        https://github.com/jamescasbon/amptools
        """


    @cmdln.option("-o", "--outfile", action="store", default='-',
                  help="Output file (default stdout)")
    @cmdln.option("-b", "--offset-bases", action="store", default=10,
                help="Allowed difference between expected start of amplicon and start of read (default 10)")
    
    def do_clip_one(self, subcmd, opts, bamfile, amplicon_location, trim_location):
        """${cmd_name}: Clip a single amplicon

            Clip reads from BAMFILE that match AMPLICON_LOCATION to TRIM_LOCATION.
            Both locations should be zero based and in the format chr:start-end. 
        
            ${cmd_usage}
            ${cmd_option_list}
        """
        
        samfile = pysam.Samfile(bamfile, 'rb')
        outfile = pysam.Samfile(opts.outfile, 'wb', template=samfile)
        stats = Stats(' '.join(sys.argv))
        
        amp_loc = Interval.from_string(amplicon_location)
        trim_loc = Interval.from_string(trim_location)
        
        if not amp_loc.contains(trim_loc):
            print('trim location not contained in amplicon location, impossible trim', file=sys.stderr)
            sys.exit(1)
        
        amplicon = Amplicon(chr=amp_loc.chrom, start=amp_loc.start, end=amp_loc.end,
            trim_start = trim_loc.start, trim_end=trim_loc.end, 
            external_id = amplicon_location, stats = stats
        )
        
        trimmed = amplicon.clipped_reads(samfile, offset_allowed=opts.offset_bases)
        map(outfile.write, trimmed)
        
        stats.report(sys.stderr)
        

    @cmdln.option("-d", "--delimiter", action="store", default=',',
                  help="Delimiter for design file (default ,)")
    @cmdln.option("-a", "--amplicon-column", action="store", default='amplicon_location',
                help="Column for amplicons in design file (default 'amplicon_location')")
    @cmdln.option("-t", "--trim-column", action="store", default='trim_location',
                help="Column for trim location in design file (default 'trim_location')")
                
    @cmdln.option("-o", "--outfile", action="store", default='-',
                  help="Output file (default stdout)")
    @cmdln.option("-b", "--offset-bases", action="store", default=10,
                help="Allowed difference between expected start of amplicon and start of read (default 10)")

    def do_clip(self, subcmd, opts, bamfile, design):
        """${cmd_name}: Clip a set of amplicons

            Clip reads from BAMFILE that match amplicon locations described in DESIGN.
            Design is expected to be delimited file that describes the amplicon and 
            trim locations.
                    
            ${cmd_usage}
            ${cmd_option_list}
        """
        stats = Stats(' '.join(sys.argv))        
        amplicons = []
        for row in csv.DictReader(file(design), delimiter=opts.delimiter): 
            amp_loc = Interval.from_string(row[opts.amplicon_column])
            trim_loc = Interval.from_string(row[opts.trim_column])

            if not amp_loc.contains(trim_loc):
                print('trim location not contained in amplicon location, impossible trim', file=sys.stderr)
                sys.exit(1)

            amplicon = Amplicon(chr=amp_loc.chrom, start=amp_loc.start, end=amp_loc.end,
                trim_start = trim_loc.start, trim_end=trim_loc.end, 
                external_id = str(amp_loc), stats = stats
            )
            amplicons.append(amplicon)
        
        samfile = pysam.Samfile(bamfile, 'rb')
        outfile = pysam.Samfile(opts.outfile, 'wb', template=samfile)
    
        for amplicon in amplicons: 
            trimmed = amplicon.clipped_reads(samfile, offset_allowed=opts.offset_bases)
            map(outfile.write, trimmed)
    
        stats.report(sys.stderr)
