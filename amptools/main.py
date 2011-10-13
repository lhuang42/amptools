from __future__ import print_function

import sys
import csv
import itertools
import cmdln
from collections import Counter, namedtuple
import pysam
from amplicon import Amplicon, load_amplicons, PILEUP_MAX_DEPTH
from ucsc import Interval    
from stats import Stats, bias_test
import vcf

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
    @cmdln.option("-i", "--id-column", action="store", default='id',
                help="Column for trim location in design file (default 'trim_location')")
    @cmdln.option("-o", "--outfile", action="store", default='-',
                  help="Output file (default stdout)")
    @cmdln.option("-b", "--offset-allowed", action="store", default=10,
                help="Allowed difference between expected start of amplicon and start of read (default 10)")

    def do_clip(self, subcmd, opts, bamfile, design):
        """${cmd_name}: Clip and filter a set of amplicons

            Clip reads from BAMFILE that match amplicon locations described in DESIGN.
            Design is expected to be delimited file that describes the amplicon and 
            trim locations.
                    
            ${cmd_usage}
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

    def do_mark(self, subcmd, opts, bamfile, design):
        """ Mark reads in bam file and optionally clip.
            
            ${cmd_usage}
            ${cmd_option_list}
        """
        samfile = pysam.Samfile(bamfile, 'rb')
        stats = Stats(' '.join(sys.argv))        
        amplicons = load_amplicons(design, stats, opts, samfile=samfile)
        outfile = pysam.Samfile(opts.outfile, 'wb', template=samfile)

        # we need to reopen the file here to get sequential access after computin the pileups
        samfile = pysam.Samfile(bamfile, 'rb')
        for read in samfile: 
            stats.reads += 1
            
            # TODO: optimisation of the list of amplicons that are considered
            for amp in amplicons:
                if amp.matches(read):
                    amp.clip(read)
                    amp.mark(read7)
            outfile.write(read) 
        
        stats.report(sys.stderr)

    
    
    def do_allele_balance(self, subcmd, opts, bamfile, position):
        """ Report allele balances at each site 
        
            ${cmd_usage}
            ${cmd_option_list}  
        """
        samfile = pysam.Samfile(bamfile, 'rb')
        position = position.split(':')
        chrom, base = position[0], int(position[1])
        
        counts = Counter()
        
        starts = {}
        ends = {}
        
        Obs = namedtuple('Observation', ['rg', 'amp', 'base'])
        
        for puc in samfile.pileup(chrom, base, base+1, max_depth=PILEUP_MAX_DEPTH):
            
            if puc.pos == base:
                for pu in puc.pileups:
                    starts[pu.alignment.qname] = pu.qpos
            
            if puc.pos == base+1:
                for pu in puc.pileups:
                    spos = starts.get(pu.alignment.qname, None)
                    if spos is None: 
                        continue                    
                    tags = dict(pu.alignment.tags)
                    insert = pu.alignment.query[spos+1:pu.qpos] or '-'
                    
                    observation = Obs(
                        rg=tags.get('RG', None), amp=tags.get('AM', None), base=insert)
                    counts[observation] += 1

            # for pu in puc.pileups:
            #     tags = dict(pu.alignment.tags)
            #     observation = tags.get('RG', None), tags.get('AM', None), pu.alignment.query[pu.qpos] 
            #     counts[observation] += 1
            
        
        alleles = sorted(set([x.base for x in counts]))
        amps = sorted(set([x.amp for x in counts]))
        rgs = sorted(set([x.rg  for x in counts]))
        
        print('rg', 'amplicon', *alleles)
        for rg in rgs: 
            for amp in amps: 
                freqs = [counts[Obs(rg=rg,amp=amp,base=b)] for b in alleles]
                print(rg, amp, *freqs)
                
        

    def do_filter_vcf(self, subcmd, opts, vcffile):
        for line in file(vcffile):
            if line.startswith('#'): 
                sys.stdout.write(line)
                
            else: 
                line = vcf.apply_filters(line)
                sys.stdout.write(line)
 
                        
                