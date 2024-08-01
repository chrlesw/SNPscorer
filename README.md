# SNPscorer

 _SNPscorer.py  is a tidied-up version of _SNPscoreFormat_v5.py ## C White  2021-2024  ##
 
Starting with a formatted list of SNP sites and a BAM file of sequences aligned to the reference genome,
     this script compiles SNP base calls for each sequence, formats them as genotype calls for each polymorphism
     and counts and outputs tables of recombinants and parentals and the accompanying statistics.

This  script has been written and tested on Apple MacOS computers. It should be relatively simple to adapt to other computers. 

Although intended as a simple tool for the analysis of amplicon sequences, this script can deal with
   genome sequencing data, with the caveat that the number of sequences & polymorphisms analysed is limited by 
   the computer's RAM. In practice this is not a constraint for analyses of Amplicon sequences, and if needed, should be 
   easy to overcome by modifying the script. 

Call as _SNPscorer.py genome1_label genome2_label polymorphisms.txt infile.sorted.bam 

#################################################################################
# Inputs:
#################################################################################
	genome1_label  (label to identify genome 1 - this is the reference genome)
	genome2_label  (label to identify genome 2)

	polymorphisms.txt 
		This is a text file with one polymorphism per line in format: CHROM POS REF ALT ; this file can easily be generated from the sorted BAM file with this command :
			callvariants.sh in=file.sorted.bam  ref=genome_ref.fa ploidy=2 minallelefraction=0.4 rarity=0.4 minreads=1000 out=varsList.vcf ; \
			     cat varsList.vcf | awk ' {print $1,$2,$4,$5}' > polymorphisms.txt
		 Note that the script only deals with single-nucleotide variants. To include >1bp deletions, the deletion will need to be split into single-bp deletions. 

callvariants.sh in=site2DR2xLer.sorted.bam  ref=~/_informatics/seqs_ref/tair10/TAIR10_chr.fa ploidy=2 minallelefraction=0.48 rarity=0.48 minreads=1000 out=varsList.vcf ; cat varsList.vcf | awk ' {print $1,$2,$4,$5}' > polys48.txt 

	infile.sorted.bam
		a sorted .bam file with the sequences mapped to the reference genome

#################################################################################
# Outputs:
#################################################################################
	<file-name>_stats.txt
		- the number of parental and recombinant sequences and the recombination rate. These are also presented taking only "clean" sequences into account (no polymorphic sites scored as "mut").
		- tables of allele counts and frequencies for each polymorphic site
	<file-name>_snpList.txt
		- one sequence per line with the nucleotide at each SNP site ("." means deletion ; "-" means non-informative (truncated sequence?)
	<file-name>_snpListFormatted.txt
		- the <file-name>_snpList.txt file with the genome1/genome2 calls for each SNP site ("mut" means neither genome1 nor genome2; "-" means non-informative (truncated sequence?))
	<file-name>_parentalSeqs.txt
		- sequences which have only genome1 alleles (no genome2 alleles) or vice-versa
	<file-name>_cleanParentalSeqs.txt
		- parental sequences filtered to remove those with one or more polymorphic sites scored as "mut" (neither genome1 nor genome2)
	<file-name>_recSeqs.txt		
		- recombinant sequences which have both genome1 and genome2 alleles
	<file-name>_cleanRecSeqs.txt
		- recombinant sequences filtered to remove those with one or more polymorphic sites scored as "mut" (neither genome1 nor genome2)

#################################################################################
# Requires:
#################################################################################

Python (tested on version 3)

samtools
	install the samtools package  (with conda...)

callvariants.sh
	This is part of the BBMap scripts package, which should be installed from sourceforge ( https://sourceforge.net/projects/bbmap/ ) and added to your $PATH

jvarkit
	see http://lindenb.github.io/jvarkit/
	A pre-compiled jar is available at https://uncloud.univ-nantes.fr/index.php/s/4sL77oWR2BFzSBH 
	install the jar file in home directory.

######################################################################################################
