"""
Amptools is a set of tools for dealing with amplicon experiments using the SAM format.
The subcommands below can be used to annotate a BAM file, find duplicates based on
molecular counters and clip primer sequences.
"""
import argparse
import sys
import annotate
import clip


parser = argparse.ArgumentParser(prog='amptools', description=sys.modules[__name__].__doc__)
parser.add_argument('--profile', action='store_true', help='run with profiling')
parser.add_argument('--verbose', '-v', action='count', help='verbosity (use -vv for debug)')
subparsers = parser.add_subparsers(help='sub-command help')

# annotate command
parser_a = subparsers.add_parser('annotate', description=annotate.annotate.__doc__,
        help='annotate a BAM file with tags')
parser_a.set_defaults(func=annotate.annotate)
parser_a.add_argument('input', type=str, help='input BAM file')
parser_a.add_argument('--output', type=str, help='output BAM file (default stdout)', default='-')

parser_a.add_argument('--adaptor', type=str, help='Adaptor in barcode/counter file.  Use B for barcode bases and M for molecular counter bases')
annotate.MidAnnotator.customize_parser(parser_a)
annotate.AmpliconAnnotator.customize_parser(parser_a)
annotate.DbrAnnotator.customize_parser(parser_a)

# duplicates command
parser_c = subparsers.add_parser('duplicates', description=annotate.duplicates.__doc__,
        help='mark duplicates based on molecular counter')
parser_c.set_defaults(func=annotate.duplicates)
parser_c.add_argument('input', type=str, help='input BAM file')
parser_c.add_argument('--output', type=str, help='output BAM file (default stdout)', default='-')

# clip command
parser_b = subparsers.add_parser('clip', description=clip.clip.__doc__,
        help='primer clip based on amplicon annotations')
parser_b.set_defaults(func=clip.clip)
parser_b.add_argument('input', type=str, help='input file')
parser_b.add_argument('--output', type=str, help='output file', default='-')

clip.AmpliconClipper.customize_parser(parser_b)

