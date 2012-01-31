#!/usr/bin env python
""" Create test data for suite

Simple script to:
    * create reference chrom
    * produce simulated reads
    * trim those reads to produce MIDs
    * run mapping to produce a BAM
"""
import os
import random
import pyfasta
import commands
import subprocess


REFLEN = 10000
REF = 'reference.fa'

MIDS = {
    'NA1': 'AAA',
    'NA2': 'CCC',
    'NA3': 'GGG',
    'NA4': 'TTT'
}

DBRS = ['AA', 'CC', 'GG', 'TT']

AMPS = [
    (100, 400, '+'),
    (200, 500, '-')
]

READS = 'reads.fa'
NREADS = 40

ADAP = "CATG%(mid)sCATG%(dbr)sCATG"

RAW_BAM = 'raw_map.bam'


def make_reference():
    """ simulate a reference sequence """
    ref = file('reference.fa', 'w')

    for x in range(1,2):
        ref.write('>chr%s\n' % x)

        for _ in range(REFLEN):
            ref.write(random.choice('ACGT'))
        ref.write('\n')



def make_reads():
    ref = pyfasta.Fasta(REF)
    out = file('reads.fa', 'w')
    rn = 0
    for amp in AMPS:

        amp_seq = ref.sequence({'chr': 'chr1', 'start': amp[0], 'stop': amp[1], 'strand': amp[2]})
        for sample, mid in MIDS.items():
            for i in range(NREADS):
                dbr = DBRS[i%len(DBRS)]
                read = ADAP % locals()
                read += amp_seq[:random.randint(len(amp_seq)/2,len(amp_seq))]
                out.write('>XXX%04d\n' %rn)
                out.write(read + '\n')
                rn += 1


def trim_reads():
    adap = ADAP % {'mid': 'NNN', 'dbr': 'NN'}
    p = subprocess.Popen(
        ['cutadapt', '-g', adap, '--wildcard-file', 'trim.txt', READS], stdout=file('trim.fa', 'w')
    )
    p.wait()

def run_map():
    p1 = subprocess.Popen('ssaha2 -output sam reference.fa trim.fa'.split(), stdout=subprocess.PIPE)
    p2 = subprocess.Popen('grep -v SSAHA'.split(), stdin=p1.stdout, stdout=subprocess.PIPE)
    p3 = subprocess.Popen('samtools view -T reference.fa -Sb -o raw_map.bam -'.split(), stdin=p2.stdout)
    p3.wait()
    os.system('samtools view -H raw_map.bam > h.txt && samtools reheader h.txt raw_map.bam > raw_map.bam_ && mv raw_map.bam_ raw_map.bam')
    os.system('rm h.txt')

if __name__ == '__main__':
    make_reference()
    make_reads()
    trim_reads()
    run_map()

#os.system('amptools annotate --mid mids.txt --trim trim.txt --amps amps.txt --dbrs trim.txt --output anno.bam_ raw_map.bam')
#os.system('samtools sort anno.bam_ anno')
#os.system('samtools index anno.bam')
#os.system('amptools clip --output clip.bam --amps amps.txt anno.bam')


