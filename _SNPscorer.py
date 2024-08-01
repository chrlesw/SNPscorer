#!/usr/bin/env   python

######################################################################################################
##     _SNPscorer.py  is a tidied-up version of _SNPscoreFormat_v5.py ## C White  2021-2024  ##
##    
##    Starting with a formatted list of SNP sites and a BAM file of sequences aligned to the reference genome,
##         this script compiles SNP base calls for each sequence, formats them as genotype calls for each polymorphism
##         and counts and outputs tables of recombinants and parentals and the accompanying statistics.
##    
##    This  script has been written and tested on Apple MacOS computers. It should be relatively simple to adapt to other computers. 
##    
##    Although intended as a simple tool for the analysis of amplicon sequences, this script can deal with
##       genome sequencing data, with the caveat that the number of sequences & polymorphisms analysed is limited by 
##       the computer's RAM. In practice this is not a constraint for analyses of Amplicon sequences, and if needed, should be 
##       easy to overcome by modifying the script. 
##    
##    Call as _SNPscorer.py genome1_label genome2_label polymorphisms.txt infile.sorted.bam 
##    
##    #################################################################################
##    # Inputs:
##    #################################################################################
##    	genome1_label  (label to identify genome 1 - this is the reference genome)
##    	genome2_label  (label to identify genome 2)
##    
##    	polymorphisms.txt 
##    		This is a text file with one polymorphism per line in format: CHROM POS REF ALT ; this file can easily be generated from the sorted BAM file with this command :
##    			callvariants.sh in=file.sorted.bam  ref=genome_ref.fa ploidy=2 minallelefraction=0.4 rarity=0.4 minreads=1000 out=varsList.vcf ; \
##    			     cat varsList.vcf | awk ' {print $1,$2,$4,$5}' > polymorphisms.txt
##    		 Note that the script only deals with single-nucleotide variants. To include >1bp deletions, the deletion will need to be split into single-bp deletions. 
##    
##    callvariants.sh in=site2DR2xLer.sorted.bam  ref=~/_informatics/seqs_ref/tair10/TAIR10_chr.fa ploidy=2 minallelefraction=0.48 rarity=0.48 minreads=1000 out=varsList.vcf ; cat ##    varsList.vcf | awk ' {print $1,$2,$4,$5}' > polys48.txt 
##    
##    	infile.sorted.bam
##    		a sorted .bam file with the sequences mapped to the reference genome
##    
##    #################################################################################
##    # Outputs:
##    #################################################################################
##    	<file-name>_stats.txt
##    		- the number of parental and recombinant sequences and the recombination rate. These are also presented taking only "clean" sequences into account (no polymorphic ##    sites scored as "mut").
##    		- tables of allele counts and frequencies for each polymorphic site
##    	<file-name>_snpList.txt
##    		- one sequence per line with the nucleotide at each SNP site ("." means deletion ; "-" means non-informative (truncated sequence?)
##    	<file-name>_snpListFormatted.txt
##    		- the <file-name>_snpList.txt file with the genome1/genome2 calls for each SNP site ("mut" means neither genome1 nor genome2; "-" means non-informative (truncated ##    sequence?))
##    	<file-name>_parentalSeqs.txt
##    		- sequences which have only genome1 alleles (no genome2 alleles) or vice-versa
##    	<file-name>_cleanParentalSeqs.txt
##    		- parental sequences filtered to remove those with one or more polymorphic sites scored as "mut" (neither genome1 nor genome2)
##    	<file-name>_recSeqs.txt		
##    		- recombinant sequences which have both genome1 and genome2 alleles
##    	<file-name>_cleanRecSeqs.txt
##    		- recombinant sequences filtered to remove those with one or more polymorphic sites scored as "mut" (neither genome1 nor genome2)
##    
##    #################################################################################
##    # Requires:
##    #################################################################################
##    
##    Python (tested on version 3)
##    
##    samtools
##    	install the samtools package  (with conda...)
##    
##    callvariants.sh
##    	This is part of the BBMap scripts package, which should be installed from sourceforge ( https://sourceforge.net/projects/bbmap/ ) and added to your $PATH
##    
##    jvarkit
##    	see http://lindenb.github.io/jvarkit/
##    	A pre-compiled jar is available at https://uncloud.univ-nantes.fr/index.php/s/4sL77oWR2BFzSBH 
##    	install the jar file in home directory.
##    
######################################################################################################


import sys, os

if (len(sys.argv) == 5) :
	genome1_label = sys.argv[1]
	genome2_label = sys.argv[2]
	polysFile = sys.argv[3]
	theBamFile = sys.argv[4]
	nameroot = sys.argv[4].split(".")[0]
	
	selSNPfileName = nameroot + "_selSnps.txt"

	print()
	print("##################################################")
	print(" _SNPscorer.py : starting with a formatted list of SNP sites")
	print("    and a sorted BAM file of sequences aligned to the reference genome,")
	print("    this script compiles SNP base calls for each sequence, formats them ")
	print("    as genotype calls, counts and outputs tables of recombinants and parentals.")
	print("")
	print(" call as _SNPscorer.py genome1_label genome2_label polymorphisms.txt infile.sorted.bam ")
	print("                                                                     C.White 2021-2024 ")
	print("##################################################")

else : 
	print("##################################################")
	print(" Exit due to missing parameter. ")
	print("  call as _SNPscorer.py genome1_label genome2_label polymorphisms.txt infile.sorted.bam ")

	print("##################################################")
	sys.exit()

########################### set up key list and polymorphism dictionary ###################################

currSeq = { "name" : '-'}
refsnps = {}

with open(polysFile) as f :
	chr = coord = val1 = val2 = key = ""
	for l in f:
		if (l[0] != "#") :
			(chr, coord, val1, val2) = l.split()
			key = chr + "_" + coord
			currSeq[key] = "-"
			refsnps[key] = val1, val2
			chr = coord = val1 = val2 = key = ""

thekeys = list(refsnps.keys())

########################################################################################################
########################### set up allele freq dictionary ###################################


numGenome1 = {}
numGenome2 = {}
numMut = {}
numNoInf = {}

for i in thekeys :
	numGenome1[i] = 0
	numGenome2[i] = 0
	numMut[i] = 0
	numNoInf[i] = 0


########################################################################################################

########################## generate input file with list of bases at the snp sites for all seqs   #################
###### this generates and the runs the command line to output the BAM file as lines of tab-separated values #######
###### and keeps only those at the polymorphism positions (taken from the input polymorphisms table)        #######

print("###  extracting all positions from the BAM file  ###")

os.system("samtools view -h " +  str(theBamFile) + "  | java -jar  ~/jvarkit.jar sam2tsv | awk ' {print $1,$6,$4,$8,$4\"_\"$8} ' > tmpOut.txt ")

print("extracting lines corresponding to the snp positions from the BAM file")
print()

theSelSnps = open(selSNPfileName, 'w')
count = 0

with open("tmpOut.txt") as f :
	chr = coord = val1 = val2 = key = ""

	for l in f: 
		if (l[0] != "#") :
			currLine = l.split()
			for currKey in thekeys : 
				if currLine[4] == currKey :
					print(currLine[0], currLine[1], currLine[2], currLine[3], currLine[4], file=theSelSnps)
					count+=1
	
print(str(count), " lines written")
theSelSnps.close()
f.close()
os.system("rm tmpOut.txt")
print()

########################################################################################################
########################## collect snp base calls per sequence  ######################################## 
print()
print("###  collecting snp base calls per sequence  ###")
print()

currSNP = currName = theName = ""
currLine = []

allSNPS = open(selSNPfileName, 'r')
firsttime = "yes"

seqsSNPs = open(nameroot + "_snpList.txt", 'w')

for theLine in allSNPS.readlines() :
	currLine = theLine.split()

	if "#" in currLine[0] :
		continue  # this skips header lines

	###  
	if currLine[0] == currName :
		currSeq["name"] = currLine[0]
		addr = currLine[4]
		currSeq[addr] = currLine[1]
	else : 
		currName = currLine[0]
		if (firsttime == "yes") :
			currSeq["name"] = currLine[0]
			addr = currLine[4]
			currSeq[addr] = currLine[1]
			for key in currSeq :
				print(key , sep='\t', end='\t',  file=seqsSNPs  )
			print( file=seqsSNPs)
			firsttime = "no"
		else:
			for key in currSeq :
				print(currSeq[key] , sep='\t', end='\t' , file=seqsSNPs)
				currSeq[key] = '-'
			print( file=seqsSNPs)
			addr = currLine[4]
			currSeq[addr] = currLine[1]

print("", file=seqsSNPs)

allSNPS.close()
seqsSNPs.close()

os.system ("rm " + selSNPfileName)

########################################################################################################
########################## now format them ############################################################# 
print()
print("###  formatting snp base calls per sequence  ###")
print()

infile = open(nameroot + "_snpList.txt", 'r')
outfile = open(nameroot + "_snpListFormatted.txt", 'w')

currSNP = currName = theName = ""
currLine = []

currSNP_numGenome1 = currSNP_numGenome2 = 0


#first print the column titles line
print("#name", end='\t' , file=outfile)
for key in refsnps :
	print(key, sep='\t', end='\t', file=outfile )

print(file=outfile)

for theLine in infile.readlines() :
	currLine = theLine.split()

	if ((len(currLine) == 0) or (currLine[0] == "name")) :
		continue #this skips directly to next readline

	count = 0
	if (len(currLine) > 1) :
		print(currLine[0], end='\t', file=outfile)

		for i in range(len(thekeys)) :
			if (i < len(currLine) - 1) :
				currSNP = currLine[i+1]  
				snpAddr = thekeys[count] 
				count +=1
				if (currSNP ==  refsnps[snpAddr][0]) : 
					print(genome1_label, end='\t' , file=outfile)
					numGenome1[snpAddr] +=1
				elif (currSNP ==  refsnps[snpAddr][1]) : 
					print(genome2_label, end='\t', file=outfile)
					numGenome2[snpAddr] +=1
				elif (currSNP ==  "-")  : #this labels as noinf those that haven't been filled in (no info in BAM)
					print("noInf", end='\t', file=outfile)
					numNoInf[snpAddr] +=1
				else  :  #this labels as mut those all others 
					print("mut", end='\t', file=outfile)
					numMut[snpAddr] +=1

		print(file=outfile)
	for i in range(len(currLine)) :
		currLine[i] = ''

print( file=outfile)

infile.close()
outfile.close()


########################################################################################################
########################### now get stats   ############################################################ 
###### (the "clean" files are filtered to exclude all sequences with one or more polymorphic sites scored as "mut")

print("###  collecting stats  ###")

infile = open(nameroot + "_snpListFormatted.txt", 'r')
rec_file = open(nameroot + "_recSeqs.txt", 'w')
parentals_file = open(nameroot + "_parentalSeqs.txt", 'w')
stats_file = open(nameroot + "_stats.txt", 'w')

clean_rec_file = open(nameroot + "_cleanRecSeqs.txt", 'w')
clean_parentals_file = open(nameroot + "_cleanParentalSeqs.txt", 'w')

currSNP = currName = theName = ""
currLine = []
r = rec = genome_1 = genome_2 = non_informative = total = 0
r_clean = clean_parentals = clean_rec  = total_clean = 0

print("#name", end='\t' , file=rec_file)
print("#name", end='\t' , file=parentals_file)
print("#name", genome1_label , genome2_label, "Rec", "r", "total", "clean_rec", "clean_parentals", "r_clean", sep='\t', file=stats_file)

print("#name", end='\t' , file=clean_rec_file)
print("#name", end='\t' , file=clean_parentals_file)

for key in refsnps :
 	print(key, sep='\t', end='\t', file=rec_file )
 	print(key, sep='\t', end='\t', file=parentals_file )

 	print(key, sep='\t', end='\t', file=clean_rec_file )
 	print(key, sep='\t', end='\t', file=clean_parentals_file )

print(file=parentals_file)
print(file=rec_file)

print(file=clean_rec_file)
print(file=clean_parentals_file)

for theLine in infile.readlines() :
	if (genome1_label in theLine):
		if (genome2_label in theLine):
			rec+=1
			print(theLine.rstrip(), file=rec_file)
			if not ( "mut" in theLine) :
				print (theLine.rstrip(), file=clean_rec_file)
				clean_rec +=1
		else: 
			genome_1 +=1
			print(theLine.rstrip(), file=parentals_file)
			if not ( "mut" in theLine) :
				print (theLine.rstrip(), file=clean_parentals_file)
				clean_parentals +=1
	elif (genome2_label in theLine):
		genome_2 +=1
		print(theLine.rstrip(), file=parentals_file)
		if not ( "mut" in theLine) :
			print (theLine.rstrip(), file=clean_parentals_file)
			clean_parentals +=1
	else :
		non_informative +=1

total = genome_1+genome_2+rec

total_clean = clean_rec + clean_parentals

if (total > 0): r = rec/total
else: r = 0

if (total_clean > 0): r_clean = clean_rec/total_clean
else: r_clean = 0


print(nameroot, str(genome_1) , str(genome_2), str(rec), str(r), str(total), str(clean_rec), str(clean_parentals), str(r_clean), sep='\t', file=stats_file)
print("###############  " + nameroot + "  ###############")
print("", genome1_label + " seqs = " + str(genome_1) , genome2_label + " seqs = " + str(genome_2), "recombinants = " + str(rec), "r = " + str(r), "total seqs = " + str(total), sep='\n')
print("Finished. Analysed " +  str(total) + " sequences")
print()
print("clean rec = " + str(clean_rec) , "clean parentals = " + str(clean_parentals) , "r_clean = " + str(r_clean) , sep='\n')
print("Analysed " +  str(total_clean) + " clean sequences (" + genome1_label + " or " + genome2_label + ", at all sites)")
print()
print(str(non_informative), "neither " + genome1_label + " nor " + genome2_label + " seqs were not scored")
print("")
print("###############################################")

#####################################################################
#############       print allele counts          ####################

print( file=stats_file )
print( file=stats_file )
print( "allele counts", file=stats_file )
print( file=stats_file )


print("", sep='\t', end='\t', file=stats_file )
for key in thekeys :
	print(key, sep='\t', end='\t', file=stats_file )
print( file=stats_file )

print("numGenome1", sep='\t', end='\t', file=stats_file )
for key in thekeys :
	print(numGenome1[key], sep='\t', end='\t', file=stats_file )
print( file=stats_file )

print("numGenome2", sep='\t', end='\t', file=stats_file )
for key in thekeys :
	print(numGenome2[key], sep='\t', end='\t', file=stats_file )
print( file=stats_file )

print("numMut", sep='\t', end='\t', file=stats_file )
for key in thekeys :
	print(numMut[key], sep='\t', end='\t', file=stats_file )
print( file=stats_file )

print("numNoInf", sep='\t', end='\t', file=stats_file )
for key in thekeys :
	print(numNoInf[key], sep='\t', end='\t', file=stats_file )
print( file=stats_file )

#####################################################################
#############       print allele frequencies     ####################

if (total > 0) :
	print( file=stats_file )
	print( file=stats_file )
	print( "allele frequencies", file=stats_file )
	print( file=stats_file )

	print("", sep='\t', end='\t', file=stats_file )
	for key in thekeys :
		print(key, sep='\t', end='\t', file=stats_file )
	print( file=stats_file )

	print("numGenome1", sep='\t', end='\t', file=stats_file )
	for key in thekeys :
		print(numGenome1[key]/total, sep='\t', end='\t', file=stats_file )
	print( file=stats_file )

	print("numGenome2", sep='\t', end='\t', file=stats_file )
	for key in thekeys :
		print(numGenome2[key]/total, sep='\t', end='\t', file=stats_file )
	print( file=stats_file )

	print("numMut", sep='\t', end='\t', file=stats_file )
	for key in thekeys :
		print(numMut[key]/total, sep='\t', end='\t', file=stats_file )
	print( file=stats_file )

	print("numNoInf", sep='\t', end='\t', file=stats_file )
	for key in thekeys :
		print(numNoInf[key]/total, sep='\t', end='\t', file=stats_file )
	print( file=stats_file )

#####################################################################
