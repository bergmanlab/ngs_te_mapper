ngs_te_mapper
=============

ngs_te_mapper.R is an R implementation of the method for detecting non-reference transposable element (TE) insertions from next-generation sequencing (NGS) data published in [Linheiro and Bergman (2012) PLoS ONE 7(2): e30008](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0030008).

Non-reference (aka de novo) TE insertions are detected using a two-stage process that relies on the presence of target site duplications (TSDs) in the region flanking the TE insertion. An overview of the two-stage mapping procedure can be found in [Figure 1](http://www.plosone.org/article/fetchObject.action?uri=info:doi/10.1371/journal.pone.0030008.g001&representation=PNG_M) of Linheiro and Bergman (2012). 

In the first stage, unaligned and unassembled whole genome shotgun sequence reads from a resequenced genome that has an integrated TE insertion are used to query against a library of reference TE sequences. Reads that span the junction of the start or end of TE and genomic flanking sequences are retained (these are often referred as 'split-reads', although in reality these reads are not split in the resequenced genome). 

In the second stage, the unique (i.e. non-TE) components of junction reads identified in the first step are aligned against a reference genome. Genomic matches are clustered and the region of overlap between sets of junction reads that span the start and end of the same reference TE are used to define the location and orientation of a non-reference TE insertions. TE insertion sites are annotated as the span of TSD and orientation is assigned in the strand field  according to a framework described in [Bergman (2012) Mob Genet Elements. 2:51-54](http://www.landesbioscience.com/journals/mge/article/19479/)

In addition to the TE insertion site mapping code, we provide an R script (ngs_te_logo.R) that clusters TSDs from the same TE family and outputs a sequence logo describing the local nucleotide preferences in a window around the TSD, as in [Figure 3](http://www.plosone.org/article/fetchObject.action?uri=info:doi/10.1371/journal.pone.0030008.g003&representation=PNG_M) of Linheiro and Bergman (2012). 

We note that the current version of ngs_te_mapper does not filter for the modal TSD length, and thus gives slightly different results to those reported in Linheiro and Bergman (2012). We also note that, as described in [Linheiro and Bergman (2012)](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0030008), the TSD logo method works best for TE families that generate fixed length TSDs (e.g. LTR retrotransposons and TIR transposons).

To run the mapping and logo methods on an example fasta or fastq file, execute:
	
	bash sourceCode/run_ngs_te_mapper.sh

For more details on how the main script (sourceCode/ngs_te_mapper.R) runs check the sourceCode/run_ngs_te_mapper.sh. 
The main script takes has input:

	genome fasta file (full path)
	TE fasta file (full path)
	fasta or fastq files of short read sequences (can either take the full path or just the file names has long has the input folder for the folders is supplied (fastaFolder)
	output folder (the folder where all the output will be written in to)

When running he main script it will look for the presence of the indexed genome and TE file in the same location has the fasta file, if not there it will create a new index in the same folder.
	
Dependencies
============

  * [R](http://cran.r-project.org/)
  * [bwa](http://bio-bwa.sourceforge.net/)

Output files and folders
============

Inside the ouput folder there will be 7 other folders:

	align_te (sam files from the alignment of the short read data to the TE file)
	aligned_te (fasta file of the selected reads)
	aligned_genome (sam files of the alignment of the selected reads to the genome file)
	bed_tsd (the files for all insertions and for the reads that predicted the new and old sites)
	metadata (the data for all insertions with the following columns - chrom;start;end;tsd_length;strand;teName;strain;nReads;insertion)
	r_data_files (contains the R workspace file from the output of ngs_te_mapper.R)
	logo (contains a pdf file of the logos for the different insertions)

And one file created by the ngs_te_logo.R script that concatenates all the .tsv files in the metadata folder
 
  

