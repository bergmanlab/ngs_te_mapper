ngs_te_mapper
=============

ngs_te_mapper is an R implementation of the method for detecting non-reference transposable element (TE) insertions from next-generation sequencing data published in [Linheiro and Bergman (2012) PLoS ONE 7(2): e30008](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0030008).

Non-reference (aka de novo) TE insertions are detected using a two-stage process that relies on the presence of target site duplications (TSDs) in the region flanking the TE insertion. In the first stage, unaligned and unassembled whole genome shotgun sequence reads from a resequenced genome that has an integrated TE insertion are used to query against a library of reference TE sequences. Reads that span the junction of the start or end of TE and genomic flanking sequences are retained (these are often referred as 'split-reads', although in reality these reads are not split in the resequenced genome). 

In the second stage, the unique (i.e. non-TE) components of junction reads identified in the first step are aligned against a reference genome. Genomic matches are clustered and the region of overlap between sets of junction reads that span the start and end of the same reference TE are used to define the location and orientation of a non-reference TE insertions.

An overview of the two-stage mapping procedure can be found in [figure 1](http://www.plosone.org/article/fetchObject.action?uri=info:doi/10.1371/journal.pone.0030008.g001&representation=PNG_M) Linheiro and Bergman (2012).

In addition to the TE insertion site mapping code, we provide an R script (ngs_te_logo.R) that clusters TSDs from the same TE family and outputs a sequence logo describing the local nucleotide preferences in a window around the TSD (as in [figure 3](http://www.plosone.org/article/fetchObject.action?uri=info:doi/10.1371/journal.pone.0030008.g003&representation=PNG_M)). 

We note that the current version of ngs_te_mapper does not filter for the modal TSD length, and thus gives slightly different results to those reported in Linheiro and Bergman (2012). We also note that, as described in [Linheiro and Bergman (2012) PLoS ONE 7(2): e30008](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0030008), the TSD logo method works best for TE families that generate fixed length TSDs (e.g. LTR retrotransposons and TIR transposons).

To run the mapping and logo methods on an example fasta file, execute:
	
	bash sourceCode/run_ngs_te_mapper.sh
	
Dependencies
============

  * [R](http://cran.r-project.org/)
  * [blat](http://hgwdev.cse.ucsc.edu/~kent/src/blatSrc.zip) 

