Introduction
------------

Next generation sequencers produce a large volume of reads which need to be
processed into relevant data about the biological samples.  There is now a set
of well supported formats and tools for storing the mapped reads (SAM) and
variant calls (VCF).  The data processing tools usually assume that the
biochemical processing of the samples is whole genome or enrichment via
hybridization.  We were interested in processing data from highly multiplexed
experiments using enrichment via PCR amplification and developed amptools to
address the particular issues we have found with this approach.

There are several points where we found the need to develop amplicon specific
tools: 

* When using overapping amplicons, primer clipping cannot happen on the raw read
  data, since synthesized DNA in one amplicon may be identical to genomic DNA in
  another.

* PCR duplicate removal cannot be based on start positions and orientation,
  since reads from the same amplicon share these properties.

* The ability to report metrics on amplicon success from the mapped reads,
  without extra metadata files

* The ability to use the distribution of variants on particular amplicons when
  considering the validity of a variant


SAM conventions  
---------------

The description of multiplexed samples in a SAM file is already supported by the
the read group (RG) tag.  We propose that we use the 'ea' tag to denote the
expected amplicon of a read.  For each expected amplicon, we add a comment line
(CO) to the sam header::

   @CO	{"ac": "chr22:42522562-42522706", "st": "-1", "type": "ea", "id": "pcr_CYP2D6_00", "tc": "chr22:42522585-42522683"}
   @CO	{"ac": "chr22:42522686-42522829", "st": "-1", "type": "ea", "id": "pcr_CYP2D6_01", "tc": "chr22:42522709-42522808"}

each of these lines is a JSON encoded description of an expected amplicon, where
we use the following tags:

* id: the id of the amplicon, matches 'ea' tags in the alignment section 

* ac: amplicon coordinates used to match the reads

* st: amplicon strand, for unidirectional reads

* tc: trim coordinates to apply when clipping

In addition, we use the alignment tag 'mc' to denote any molecular counter
observed with the reads.  We use the usual tag 'BC' to store any multiplexing
barcode information.

Note, you might wonder why we are using user (lowercase) tags, but at the
moment the attempt to standardise tags was not supported on the SAM development
list.  Of course, if you wanted to add your support there...

http://sourceforge.net/mailarchive/forum.php?thread_name=CAFy2%2ByuFnri%3DwVz-XKEBZtPecYY6OaaeZLnRX6DuJPZaPpPyfw%40mail.gmail.com&forum_name=samtools-devel
    

Annotation of SAM files 
-----------------------

amptools offers a way of annotating sam files using these tags by using the
`annotate` command.  This command offers several different annotators that can
be activated and can operate on a file or a stream.  Input can be unsorted :: 

    usage: amptools annotate [-h] [--output OUTPUT] [--adaptor ADAPTOR]
                             [--rgs RGS] [--bcs-read BCS_READ]
                             [--rgs-read RGS_READ] [--library LIBRARY]
                             [--platform PLATFORM] [--exclude-rg EXCLUDE_RG]
                             [--offbyone] [--amps AMPS] [--id-column ID_COLUMN]
                             [--amplicon-column AMPLICON_COLUMN]
                             [--trim-column TRIM_COLUMN] [--delimiter DELIMITER]
                             [--offset-allowed OFFSET_ALLOWED]
                             [--exclude-offtarget] [--clip] [--counters COUNTERS]
                             input

    Annotate reads in a SAM file with tags. Use one or more available annotators
    below to add tags to a SAM file.

    positional arguments:
      input                 input BAM file

    optional arguments:
      -h, --help            show this help message and exit
      --output OUTPUT       output BAM file (default stdout)
      --adaptor ADAPTOR     Adaptor in barcode/counter file. Use B for barcode
                            bases and M for molecular counter bases

    RG annotation:
      Annotate BAM file with read groups (RGs) based on molecular barcodes. This
      annotator adds RGs to the header and assigns each read a RG tag with the
      reag group and a BC tag with the barcode read.

      --rgs RGS             file containing whitespace separated BC, RG pairs
                            (enables annotator)
      --bcs-read BCS_READ   file containing whitespace separated BC, read
                            accession pairs
      --rgs-read RGS_READ   file containing whitespace separated RG, read
                            accession pairs
      --library LIBRARY     (optional) library to use in RG header
      --platform PLATFORM   (optional) platform to use in RG header
      --exclude-rg EXCLUDE_RG
                            (optional) platform to use in RG header
      --offbyone            Allow off by one errors in the MID

    EA annotation:
      Annotate reads that match expected amplicons (EAs). Mark each read that
      matches and expected amplicon with an EA tag and add a EA section to the
      header.

      --amps AMPS           Delimited file describing amplicons (required)
      --id-column ID_COLUMN
                            ID column (default id)
      --amplicon-column AMPLICON_COLUMN
                            Amplicon coordinate column (default amplicon)
      --trim-column TRIM_COLUMN
                            Amplicon trim coordinates (default trim)
      --delimiter DELIMITER
                            file delimiter (default TAB)
      --offset-allowed OFFSET_ALLOWED
                            Allowed bases between read start and amplicon start
                            (default 10)
      --exclude-offtarget   Exclude reads not matching target amplicons
      --clip

    MC annotation:
      Annotate BAM file with molecular counters (MCs). This annotator adds a MC
      tag for each read contaning any molecular counter sequence read.

      --counters COUNTERS   File containing whitespace separated MC, read accesion


Read group annotation 
.....................

This annotator adds read groups to the header and alignments.  It can optionally
add barcodes to the reads, if provided.  You can activate it with `--rgs` which
takes a file containing a list of barcodes and read groups.  You must then
specify either `--rgs-read` containing the read group and accession for each
read or `--bcs-read` which contains the barcode read and accession.  When
providing the barcodes read the default strategy is to expect exact matching
barcodes.  You can change the barcode matching strategy by using `--offbyone`
which precomputes off by one errors for all barcodes or `--ngram` which uses an
ngram score to find the closest matching barcode.  You can add extra metadata
to the RG header lines using the `--library` and `--platform` flags.

Expected amplicon annotation
............................

The amplicon annotation is activated by the `--amps` flag which requires a
delimited file that contains the ids, amplicon coordinates and trim coordinates.
You can use the `--offset-allowed` flag to control the mismatch allowed between
expected and observed start positions.  Use `--exclude-offtarget` to remove
reads not matching expected amplicons.  Use `--clip` to clip to the trim
coordinates.  Short reads that only contain primer sequence will be excluded
from the output.

Molecular Counters 
..................

Molecular counters can be added with the `--




