from __future__ import print_function
import argparse

import annotate
import clip


parser = argparse.ArgumentParser(prog='amptools')
subparsers = parser.add_subparsers(help='sub-command help')


parser_a = subparsers.add_parser('annotate', help='annotate a bam file')
parser_a.add_argument('input', type=str, help='input file')
parser_a.add_argument('--output', type=str, help='output file', default='-')
annotate.MidAnnotator.customize_parser(parser_a)
annotate.AmpliconAnnotator.customize_parser(parser_a)
annotate.DbrAnnotator.customize_parser(parser_a)
parser_a.set_defaults(func=annotate.annotate)

parser_b = subparsers.add_parser('clip', help='clip help')
parser_b.add_argument('input', type=str, help='input file')
parser_b.add_argument('--output', type=str, help='output file', default='-')
parser_b.set_defaults(func=clip.clip)
clip.AmpliconClipper.customize_parser(parser_b)

