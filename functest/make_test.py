import os
import random
import pyfasta
import commands
import subprocess


REFLEN = 1000
REF = 'reference.fa'

MIDS = {
    'NA1': 'AAA',
    'NA2': 'CCC',
    'NA3': 'GGG',
    'NA4': 'TTT'
}

AMPS = [
    (100, 400, '+'),
    (200, 500, '-')
]

READS = 'reads.fa'
NREADS = 40

ADAP = "CATG%(mid)sCATG"


def make_reference():
    """ simulate a reference sequence """
    ref = file('reference.fa', 'w')

    ref.write('>chr1\n')

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
            for _ in range(NREADS):
                read = ADAP % locals()
                read += amp_seq[:random.randint(len(amp_seq)/2,len(amp_seq))]
                out.write('>XXX%04d\n' %rn)
                out.write(read + '\n')
                rn += 1


def trim_reads():
    adap = ADAP % {'mid': 'NNN'}
    p = subprocess.Popen(
        ['cutadapt', '-g', adap, '--wildcard-file', 'trim.txt', READS], stdout=file('trim.fa', 'w')
    )
    p.wait()

def run_map():
    p1 = subprocess.Popen('ssaha2 -output sam reference.fa trim.fa'.split(), stdout=subprocess.PIPE)
    p2 = subprocess.Popen('grep -v SSAHA'.split(), stdin=p1.stdout, stdout=subprocess.PIPE)
    p3 = subprocess.Popen('samtools view -T reference.fa -Sb -o raw_map.bam -'.split(), stdin=p2.stdout)
    p3.wait()


exists = os.path.exists
if not exists(REF): make_reference()
if not exists(READS): make_reads()
trim_reads()
run_map()
