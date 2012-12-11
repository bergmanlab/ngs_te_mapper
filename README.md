ngs_te_mapper
=============

ngs_te_mapper is an R implementation of the method for detecting non-reference transposable element (TE) insertions from next-generation sequencing data published in [Linheiro and Bergman (2012) PLoS ONE 7(2): e30008](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0030008)

Non-reference (aka de novo) TE insertions are detected using a two-stage process that relies on the presence of target site duplications (TSDs) in the region flanking the TE insertion. In the first stage, unaligned and unassembled whole genome shotgun sequence reads from a resequenced genome that has an integrated TE insertion are used to query against a library of reference TE sequences. Reads that span the junction of the start or end of TE and genomic flanking sequences are retained (these are often referred as 'split-reads', although in reality these reads are not split in the resequenced genome). 

In the second stage, the unique (i.e. non-TE) components of junction reads identified in the first step are aligned against a reference genome. Genomic matchs are clustered and the region of overlap between sets of junction reads that span the start and end of the same reference TE are used to define the location and orientation of a non-reference TE insertions.

Dependencies
============

  * [R](http://cran.r-project.org/)
  * [blat](http://hgwdev.cse.ucsc.edu/~kent/src/blatSrc.zip) 
  * [seqtk](https://github.com/lh3/seqtk)

