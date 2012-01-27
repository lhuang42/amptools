import os
import unittest
import pysam
import tempfile
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
        args.mids = path_to('mids.txt')
        args.mids_read = path_to('trim.txt')
        args.platform = 'LS454'
        args.library = 'L1'

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
            assert 'RG' in dict(r.tags)
            assert 'BC' in dict(r.tags)


    def test_DbrAnnotator(self):
        sf = raw_bam()
        header = sf.header
        args = MockArgs()
        args.counters = path_to('trim.txt')

        anno = annotate.DbrAnnotator(args, header)

        # each read should have a RG
        for r in sf:
            anno(r)
            assert 'MC' in dict(r.tags)


    def test_AmpliconAnnotator(self):
        sf = raw_bam()
        header = sf.header
        args = MockArgs()
        args.amps = path_to('amps.txt')
        args.delimiter = '\t'
        args.id_column = 'id'
        args.amplicon_column = 'amplicon'
        args.trim_column = 'trim'
        args.offset_allowed = 10
        args.clip = False

        anno = annotate.AmpliconAnnotator(args, header)

        # the MID annotator should put the RGs into the header
        print sf.header
        self.assertEquals('EA' in header, True)
        self.assertEquals(len(header['EA']), len(make_test.AMPS))
        for d in header['EA']:
            print d
            assert d['ID'] in 'AB'
            assert d['AC']
            assert d['TC']
            assert d['ST']

        # each read should have an EA
        for r in sf:
            anno(r)
            assert 'EA' in dict(r.tags)

    def test_duplicates(self):
        tmp = tempfile.mktemp()
        tmpo = tempfile.mktemp()

        print 'using tempfile', tmp
        try:
            # create dbr marked file
            os.system('amptools annotate --output %s --mids %s --mids-read %s --counters %s %s' % (
                tmp, path_to('mids.txt'), path_to('trim.txt'), path_to('trim.txt'), path_to(make_test.RAW_BAM)))

            args = MockArgs()
            args.input = tmp
            args.output = tmpo

            annotate.duplicates(args)

            non_dups = filter(lambda x: not x.is_duplicate, pysam.Samfile(tmpo))
            assert len(non_dups) == len(make_test.MIDS) * len(make_test.AMPS)

        finally:
            os.unlink(tmp)
            os.unlink(tmpo)


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
                a = dict(r.tags)['EA']
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

