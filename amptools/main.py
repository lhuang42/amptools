from __future__ import print_function

import sys
import csv
import itertools
import cmdln
from collections import Counter, namedtuple
import pysam
from amplicon import Amplicon, load_amplicons, PILEUP_MAX_DEPTH
from ucsc import Interval    
from stats import Stats


class Amptools(cmdln.Cmdln):
    """Usage:
            amptools SUBCOMMAND [ARGS...]
            amptools help SUBCOMMAND

        Amptools provides a set of tools for handling amplicon experiments.

        ${command_list}
        ${help_list}

        For additional information, bug reports and development, see 
        https://github.com/jamescasbon/amptools
        """

    @cmdln.option("-d", "--delimiter", action="store", default=',',
                  help="Delimiter for design file (default ,)")
    @cmdln.option("-a", "--amplicon-column", action="store", default='amplicon_location',
                help="Column for amplicons in design file (default 'amplicon_location')")
    @cmdln.option("-t", "--trim-column", action="store", default='trim_location',
                help="Column for trim location in design file (default 'trim_location')")
    @cmdln.option("-i", "--id-column", action="store", default='id',
                help="Column for amplicon identifier in design file (default 'external_id')")
    @cmdln.option("-o", "--outfile", action="store", default='-',
                  help="Output file (default stdout)")
    @cmdln.option("-b", "--offset-allowed", action="store", default=10,
                help="Allowed difference between expected start of amplicon and start of read (default 10)")

    def do_clip(self, subcmd, opts, bamfile, amplicons):
        """${cmd_name}: Find and clip reads matching amplicons.
        
            Find reads from amplicons in the input and write clipped reads to 
            output.  Writes the AM tag for matches.  Use 'mark' if you want all 
            input reads (including non matching) in the output.        
        
            ${cmd_usage}
            BAMFILE: input reads (use - for stdin)
            AMPLICONS: a file listing amplicons and trim locations.

            ${cmd_option_list}
        """
        stats = Stats(' '.join(sys.argv))  
        opts.clip = False      
        samfile = pysam.Samfile(bamfile, 'rb')
        amplicons = load_amplicons(design, stats, opts)        
        outfile = pysam.Samfile(opts.outfile, 'wb', template=samfile)
    
        for amplicon in amplicons: 
            trimmed = amplicon.clipped_reads(samfile, mark=True)
            map(outfile.write, trimmed)
    
        stats.report(sys.stderr)
        
    @cmdln.option("-d", "--delimiter", action="store", default=',',
                  help="Delimiter for design file (default ,)")
    @cmdln.option("-a", "--amplicon-column", action="store", default='amplicon_location',
                help="Column for amplicons in design file (default 'amplicon_location')")
    @cmdln.option("-t", "--trim-column", action="store", default='trim_location',
                help="Column for trim location in design file (default 'trim_location')")
    @cmdln.option("-i", "--id-column", action="store", default='id',
                help="Column for trim location in design file (default 'trim_location')")
    @cmdln.option("-o", "--outfile", action="store", default='-',
                  help="Output file (default stdout)")
    @cmdln.option("-b", "--offset-allowed", action="store", default=10,
                help="Allowed difference between expected start of amplicon and start of read (default 10)")
    @cmdln.option("-c", "--clip", action="store_true", 
                help="Trim reads that match amplicons")
                
    def do_mark(self, subcmd, opts, bamfile, amplicons):
        """${cmd_name}: Mark reads matching amplicons and optionally clip.
            
            Walk a BAM file and mark any matching amplicons using the AM tag.
            Outputs a modified BAM.  Use 'clip' if you want only reads matching 
            amplicons in the output.
            
            ${cmd_usage}
            BAMFILE: input reads (use - for stdin)
            AMPLICONS: a file listing amplicons and trim locations.
            
            ${cmd_option_list}
        """
        samfile = pysam.Samfile(bamfile, 'rb')
        stats = Stats(' '.join(sys.argv))        
        amplicons = load_amplicons(design, stats, opts, samfile=samfile)
        outfile = pysam.Samfile(opts.outfile, 'wb', template=samfile)

        # we need to reopen the file here to get sequential access after computin the pileups
        samfile = pysam.Samfile(bamfile, 'rb')
        for read in samfile: 
            
            # TODO: optimisation of the list of amplicons that are considered
            for amp in amplicons:
                if amp.matches(read):
                    amp.clip(read)
                    amp.mark(read)
            outfile.write(read) 
        
        stats.report(sys.stderr)

    
        
