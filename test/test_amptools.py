import os
import unittest
import pysam
import json
import tempfile
import commands
from collections import Counter

import make_test
from amptools import annotate
from amptools import clip

def path_to(testfile):
    op = os.path
    return op.join(op.dirname(__file__), testfile)

def raw_bam():
    return pysam.Samfile(path_to(make_test.RAW_BAM))


class MockArgs(object):
    pass


class AnnotateTest(unittest.TestCase):

    def test_MidAnnotator(self):

        sf = raw_bam()
        header = sf.header
        args = MockArgs()
        args.rgs = path_to('mids.txt')
        args.bcs_read = path_to('trim.txt')

        args.platform = 'LS454'
        args.library = 'L1'
        args.exclude_rg = None
        args.offbyone = None
        args.ngram = None

        anno = annotate.MidAnnotator(args, header)

        # the MID annotator should put the RGs into the header
        print sf.header
        self.assertEquals('RG' in header, True)
        self.assertEquals(len(header['RG']), len(make_test.MIDS))
        for d in header['RG']:
            print d
            assert d['ID'] in make_test.MIDS
            assert d['PL'] == 'LS454'
            assert d['LB'] == 'L1'

        # each read should have a RG
        for r in sf:
            anno(r)
            print dict(r.tags)
            assert 'RG' in dict(r.tags)
            assert 'BC' in dict(r.tags)


    def test_DbrAnnotator(self):
        sf = raw_bam()
        header = sf.header
        args = MockArgs()
        args.counters = path_to('trim2.txt')
        anno = annotate.DbrAnnotator(args, header)

        # each read should have a RG
        for r in sf:
            anno(r)
            assert annotate.TAG_COUNT in dict(r.tags)
            print dict(r.tags), make_test.DBRS
            assert dict(r.tags)[annotate.TAG_COUNT] in make_test.DBRS


    def test_AmpliconAnnotator(self):
        sf = raw_bam()
        header = sf.header
        args = MockArgs()
        args.input = sf
        args.amps = path_to('amps.txt')
        args.delimiter = '\t'
        args.id_column = 'id'
        args.amplicon_column = 'amplicon'
        args.trim_column = 'trim'
        args.offset_allowed = 10
        args.clip = False
        args.exclude_offtarget = False

        anno = annotate.AmpliconAnnotator(args, header)

        # the MID annotator should put the RGs into the header
        print sf.header
        self.assertEquals('CO' in header, True)
        self.assertEquals(len(header['CO']), len(make_test.AMPS))

        for d in header['CO']:
            d = json.loads(d)

            assert d['id'] in 'AB'
            assert d['ac']
            assert d['tc']
            assert d['st']

        # each read should have an EA
        for r in sf:
            anno(r)
            assert annotate.TAG_AMP in dict(r.tags)

    def test_AmpliconAnnotator_with_clip(self):
        sf = raw_bam()
        header = sf.header
        args = MockArgs()
        args.amps = path_to('amps.txt')
        args.delimiter = '\t'
        args.id_column = 'id'
        args.amplicon_column = 'amplicon'
        args.trim_column = 'trim'
        args.offset_allowed = 10
        args.clip = True
        args.input = sf
        args.exclude_offtarget = False

        anno = annotate.AmpliconAnnotator(args, header)

        # each read should have an EA
        for r in sf:
            anno(r)
            a = dict(r.tags).get(annotate.TAG_AMP, None)

            assert a is not None

            start, end = r.positions[0][1], r.positions[-1][1]
            if a == 'A':
                assert start == 120
                assert end <= 380

            elif a == 'B':
                assert start >= 220
                assert end == 480 - 1 #TODO: check off by one





    def test_duplicates(self):
        tmp = tempfile.mktemp()
        tmpo = tempfile.mktemp()

        print 'using tempfile', tmp, tmpo
        try:
            # create dbr marked file
            os.system('amptools annotate --output %s --rgs %s --bcs-read %s --counters %s %s' % (
                tmp, path_to('mids.txt'), path_to('trim.txt'), path_to('trim2.txt'), path_to(make_test.RAW_BAM)))

            args = MockArgs()
            args.input = tmp
            args.output = tmpo

            annotate.duplicates(args)

            non_dups = filter(lambda x: not x.is_duplicate, pysam.Samfile(tmpo))
            print len(non_dups)
            assert len(non_dups) == len(make_test.MIDS) * len(make_test.AMPS) * len(make_test.DBRS)

        finally:
            pass
            #os.unlink(tmp)
            #os.unlink(tmpo)


class ClipTest(unittest.TestCase):

    def test_clip(self):
        tmp = tempfile.mktemp()
        tmpo = tempfile.mktemp()

        print 'using tempfile', tmp
        try:
            # create dbr marked file
            os.system('amptools annotate --output %s --amps %s %s' % (
                tmp, path_to('amps.txt'), path_to(make_test.RAW_BAM)))

            os.system('samtools sort %s %s.sort' %(tmp, tmp))
            os.system('samtools index %s.sort.bam' % (tmp))

            args = MockArgs()
            args.input = tmp + '.sort.bam'
            args.output = tmpo
            clip.clip(args)

            n = Counter()
            for r in pysam.Samfile(tmpo):
                a = dict(r.tags)[annotate.TAG_AMP]
                n[a] += 1

                start, end = r.positions[0][1], r.positions[-1][1]
                if a == 'A':
                    assert start == 120
                    assert end <= 380

                elif a == 'B':
                    assert start >= 220
                    assert end == 480 - 1 #TODO: check off by one

        finally:
            pass


expected_stats = """total 320 reads, on target 320, uniq 32
on target 100.00%
on target reads per counter: 10.00
rg,lib,amp,unique,reads
NA1,None,A,4,40
NA1,None,B,4,40
NA2,None,A,4,40
NA2,None,B,4,40
NA3,None,A,4,40
NA3,None,B,4,40
NA4,None,A,4,40
NA4,None,B,4,40"""

class StatsTest(unittest.TestCase):
    def test_stats(self):
        tmp1, tmp2 =  tempfile.mktemp(),  tempfile.mktemp()

        os.system('amptools annotate --output %s --amps %s --rgs %s --bcs-read %s --counters %s %s' % (
                tmp1, path_to('amps.txt'), path_to('mids.txt'), path_to('trim.txt'), path_to('trim2.txt'), path_to(make_test.RAW_BAM)))


        os.system('samtools sort %s %s.sort' %(tmp1, tmp1))
        os.system('samtools index %s.sort.bam' % (tmp1))

        os.system('amptools duplicates --output %s %s.sort.bam' % (tmp2, tmp1))
        print tmp1, tmp2

        statsout = commands.getoutput('amptools coverage %s' % tmp2)
        statsout = statsout.replace('\r', '')
        assert statsout == expected_stats



