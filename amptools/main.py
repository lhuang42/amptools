from __future__ import print_function
import argparse

import annotate

parser = argparse.ArgumentParser(prog='amptools')
subparsers = parser.add_subparsers(help='sub-command help')


parser_a = subparsers.add_parser('annotate', help='annotate a bam file')
parser_a.add_argument('input', type=str, help='input file')
annotate.MidAnnotator.customize_parser(parser_a)
annotate.AmpliconAnnotator.customize_parser(parser_a)
parser_a.set_defaults(func=annotate.annotate)



parser_b = subparsers.add_parser('b', help='b help')
parser_b.add_argument('--baz', choices='XYZ', help='baz help')

