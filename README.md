ngs_te_mapper
=============

ngs_te_mapper is an R implementation of the method for detecting non-reference transposable element (TE) insertions from next-generation sequencing (NGS) data published in [Linheiro and Bergman (2012) PLoS ONE 7(2): e30008](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0030008). 

Non-reference (aka _de novo_) TE insertions are detected using a two-stage process that relies on the presence of target site duplications (TSDs) in the region flanking the TE insertion. An overview of the two-stage mapping procedure can be found in [Figure 1](http://www.plosone.org/article/fetchObject.action?uri=info:doi/10.1371/journal.pone.0030008.g001&representation=PNG_M) of Linheiro and Bergman (2012). 

In the first stage, raw reads from a whole genome shotgun sequence are used to query against a library of reference TE sequences. 'Junction reads' that span the start/end of TE and genomic flanking sequences are retained. Such reads are often referred as 'split-reads', although in reality these reads are not split in the resequenced genome. 

In the second stage, the unique (i.e. non-TE) components of junction reads are aligned to a reference genome. Junction reads in the genome are clustered according to whether the unique regions are 5' or 3' to a putative TE insertion site in the genome. The region of overlap between 5' and 3' clusters of junction reads defines the location of the target site duplication (TSD) for non-reference TE insertions. 

Non-reference TE insertion sites are annotated as the span of TSD on zero-based, half-open coordinates and orientation is assigned in the strand field, following the framework described in [Bergman (2012) Mob Genet Elements. 2:51-54](http://www.landesbioscience.com/journals/mge/article/19479/). The relative orientation of the TE in junction reads is used to annotate the strand of non-reference TE insertions. 

In addition to the TE insertion site mapping code, we provide an R script (ngs_te_logo.R) that combines TSD information from insertions of the same TE family and outputs a sequence logo describing the local nucleotide preferences in a window around the TSD, as in [Figure 3](http://www.plosone.org/article/fetchObject.action?uri=info:doi/10.1371/journal.pone.0030008.g003&representation=PNG_M) of Linheiro and Bergman (2012). As described in the paper, the TSD logo method only works for TE families that generate a fixed length TSD (e.g. LTR retrotransposons and TIR transposons). 

We note that the current version of ngs_te_logo.R does not currently filter for the modal TSD length, and thus gives slightly different results to those reported in Linheiro and Bergman (2012). Additionally, the current implementation uses [BWA](http://bio-bwa.sourceforge.net/) as a short read mapping engine instead of BLAT, as was used originally in Linheiro and Bergman (2012). 


Example
=======

To run the mapping and logo methods on an test dataset, clone this repository and execute the test script as follows:

```
git clone https://github.com/bergmanlab/ngs_te_mapper.git
cd ngs_te_mapper
bash sourceCode/run_ngs_te_mapper.sh
```

This test script runs the main script `ngs_te_mapper.R` which takes five required arguments as input:
- sourceCodeFolder (full path to location of `ngs_te_mapper.R` script)
- genome (full path to reference genome fasta file)
- teFile (full path TE fasta file )
- sample (full path to fasta or fastq files of short read sequences; alternatively list of file names if the input folder option is specified, see below)
- output (full path to output folder )

`ngs_te_mapper.R` takes two additional options arguments as input:
- mapq (filter , default=60)
- tsd ( , default=20)


When running the main script it will look for the presence of the indexed genome and TE file in the same location as the genome fasta file, if not there it will create a new genome index in the same folder.
	
Output files and folders
============

ngs_te_mapper creates a main output directory called 'analysis'. Inside this directory there will be 7 other directories:
- aligned_te (.sam files of short reads aligned to the TE file, one per input sample)
- selected_reads (.fasta file of the selected reads, one file containing selected reads from all samples)
- aligned_genome (.sam file of selected reads aligned to reference genome, one file containing selected reads from all samples)
- bed_tsd (.bed file, one file containing predicted reference and non-reference TE insertions)
- logo (.pdf file, one file with one per logo per family derived from non-reference insertions)

Dependencies
============

  * [R](http://cran.r-project.org/)
  * [BWA](http://bio-bwa.sourceforge.net/)
