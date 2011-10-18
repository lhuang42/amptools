Amptools
========

Amptools is a set of tools for describing and analysing amplicon experiments 
sequenced using massively parallel DNA sequencers.  The tools are generally 
based around the [SAM][1] format and layer in some extra information that 
reflects the needs of an amplicon sequencing experiment.

Amplicon experiments
--------------------

Amplicon sequencing experiments are based around amplification of specific 
genomic regions which are then pooled for sequencing.  A given experiment 
may involve multiple amplicons and multiple samples.  We want to be able 
to answer questions such as 

* What is the distribution of the number of reads per amplicon?
* Are any observed variants seen on multiple amplicons, or are particular amplicons biased?

Amptools can help with these questions and typical steps involved in 
processing an amplicon experiment.


Primer trimming
---------------

A simple strategy for removing primers from sequencing reads looks for 
reads with the primers and removes the primer sequence.  However, in a tiled
amplicon experiment, primer sequence for one amplicon may be genomic sequence
for another amplicon.  Amptools therefore offers a method for primer trimming
mapped reads, where the actual amplicon that produced the read can be 
identified unambiguously.


PCR duplicate detection
-----------------------

PCR duplicates are usually detected in sequencing experiments by finding reads
with the same start position and orientation. Since all reads from the same 
amplicon match this criteria another approach must be used.  Amptools supports 
the inclusion of [random degenerate bases in the adaptor][2] which can be 
used to detect PCR duplicates.

SAM format tags
---------------

Amptools uses two [SAM][1] read tags to add information to a SAM file:

* The amplicon name are stored in *AM* 
* Degenerate bases are stored in *DB* 


[1]: http://samtools.sourceforge.net/SAM1.pdf
[2]: http://nar.oxfordjournals.org/content/early/2011/04/13/nar.gkr217.abstract

