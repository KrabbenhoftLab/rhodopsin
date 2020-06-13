#!/bin/sh
#Script for running genotyping analyses on nanopore Flongle data
#KM Eaton, May 2020, University at Buffalo

#ABOUT THIS SCRIPT:
#This script takes raw fastq files of barcoded amplicons from an Oxford Nanopore sequencing run. This script assumes that basecalling has already been done on these files.
#It demultiplexes reads (identifies those specific to each barcode), maps them to a reference, and calls SNPs.
#This script assumes you have the following programs installed: 
#Guppy V3.2.4 or later (available from Oxford Nanopore Technologies) 
#bwa V0.7.17 or later (Li 2013, available from https://sourceforge.net/projects/bio-bwa/files/) 
#samtools V1.9 or later (Li et al. 2009, available from https://www.htslib.org/download/) 
#bcftools V1.9 or later (Li 2010, Li 2011, Danecek et al. 2014, available from https://www.htslib.org/download/) 
#VCFtools V0.1.17 or later (Danecek et al. 2011, available from https://vcftools.github.io/downloads.html)

#USER INPUT: Keep the terms on the left, but fill in the quotes on the right according to the directions as commented below. They're filled in as an example right now.
guppydir="/home/krablab/Documents/ont-guppy-cpu/bin" #Directory containing guppy_barcoder. Make sure there's NO forward slash at the end of the path.
fastq_pass="/mnt/sdb3/Katie_Cisco/test-bash-2/fastq_pass" #Path to directory containing fastq files that have passed quality score checks (done automatically by MinKNOW/Guppy, will be located in folder called fastq_pass). Make sure there's NO forward slash at the end of the path.
barcoded_location="/mnt/sdb3/Katie_Cisco/test-bash-2/barcoded-reads" #New directory to which trimmed and demultiplexed reads should be saved. Make sure there's NO forward slash at the end of the path.
barcode_kit="EXP-PBC001" #Identity of Oxford Nanopore barcoding kit used to assign sample-specific barcodes
barcode_number="01 02 03 04 05 06 07 08 09 10 11 12" #Numbers of each barcode used in your sequencing run. We used ONT barcodes 01-12, but you may have up to 96.
db_prefix="Clavaretus_genome" #Name of database to be created by bwa. I usually just use the name of the reference genome, wihtout the file extension.
reference_genome_path="/mnt/sdb3/Katie_Cisco/test-bash-2" #Path to directory containing desired reference genome. Make sure there's NO forward slash at the end of the path.
ref_genome="Clavaretus_genome.fasta" #Filename of reference genome to which reads will be mapped. Must be in fasta format. 
all_samples="/mnt/sdb3/Katie_Cisco/test-bash-2/test-all-samples.txt" #Path to text file containing the names of each of the sorted bam files output from STEP 4. Each filename should be on a new line. 
#For example, your all_samples file should look like the following (without the hashtags indicating the comment)
#barcode01_sorted_reads.bam
#barcode02_sorted_reads.bam
#barcode03_sorted_reads.bam
#barcode04_sorted_reads.bam
#barcode05_sorted_reads.bam
#barcode06_sorted_reads.bam
#And so on, until you've included all the barcodes you've used.
pop_1="/mnt/sdb3/Katie_Cisco/test-bash-2/test-pop-1.txt" #Path to text file similar to the all_samples file described above. This file should contain the same information, but only for samples belonging to population 1. 
#For example, if I had sequenced individuals with barcodes 01 through 06, and individuals barcoded with 01, 03, and 04 belonged to population 1, the file would look like:
#barcode01_sorted_reads.bam
#barcode03_sorted_reads.bam
#barcode04_sorted_reads.bam
pop_2="/mnt/sdb3/Katie_Cisco/test-bash-2/test-pop-2.txt" #Path to text file in the same format as pop_1, but now with only info for samples belonging to population 2. 
#For example, following for as described above, if individuals with barcodes 02, 05, and 06 belong to population 2, the file would look like:
#barcode02_sorted_reads.bam
#barcode05_sorted_reads.bam
#barcode06_sorted_reads.bam

#STEP 1: Sort raw reads from fastq files by barcode, then trim barcodes and adapters
${guppydir}/guppy_barcoder --input_path ${fastq_pass} --save_path ${barcoded_location} --barcode_kits ${barcode_kit} --trim_barcodes
#Output will be X+1 number of folders in barcoded_location, where X is the number of barcodes in the kit you used.
#Each folder will be labelled with a barcode identity (barcode01, barcode02, etc.). These correspond to the barcodes you used. 
#In each folder, there will be one or more fastq files that contain reads for which that specific barcode was detected. The --trim_barcodes option, if invoked, will remove the barcode, making downstream analyses straightforward.
#There will also be a folder labelled "unclassified", which contains fastq files for reads that couldn't be assigned to a specific barcode. You can ignore this folder.

#STEP 2: Merge all your fastq files for a single barcode, rename this file with the barcode ID, and save it in a new directory within the barcoded_location folder.
cd ${barcoded_location}
mkdir ${barcoded_location}/fastq_summary/ #This creates the directory that your output files will be sent to.
for bc in ${barcode_number};
do
	cd ${barcoded_location}/barcode${bc}/
	cat *.fastq > barcode${bc}.fastq
	mv barcode${bc}.fastq ../fastq_summary/
done
#Output will be X fastq files, where X is the number of barcodes listed in barcode_number. These fastq files will be saved to a new folder in the barcoded_location directory, under the new directory fastq_summary.

#STEP 3: Map reads specific to each barcode to a reference genome using bwa.
cd ${barcoded_location}
bwa index -a bwtsw -p ${db_prefix} ${reference_genome_path}/${ref_genome} #This indexes your reference genome.
mkdir ${barcoded_location}/mapped_reads/ #This creates the directory that your output files will be sent to.
mv ${reference_genome_path}/${ref_genome} ${barcoded_location}/fastq_summary/	 #This moves your reference genome to the right directory.
mv ${barcoded_location}/${db_prefix}* ${barcoded_location}/fastq_summary/ #This moves the databases created from the reference genome to the right directory.
cd ${barcoded_location}/fastq_summary/
for bc in ${barcode_number};
do
	bwa mem -t 1 ${db_prefix} barcode${bc}.fastq > barcode${bc}_mapped_reads.sam 
	mv barcode${bc}_mapped_reads.sam ../mapped_reads/
done
#Output will be X sam files, where X is the number of barcodes listed in barcode_number. These files will be saved to the same location as in the previous step - the barcoded_location/fastq_summary directory.

#STEP 4: Use samtools to convert between filetypes, and sort and index your mapped reads.
#The -@ option invoked below is optional, it simply specifies the number of threads. All the other options are required.
mkdir ${barcoded_location}/sorted_bam_files/ #This creates the directory that your output files will be sent to.
cd ${barcoded_location}/mapped_reads/
for bc in ${barcode_number};
do
	samtools view -b -@ 1 barcode${bc}_mapped_reads.sam | samtools sort -@ 1 -o barcode${bc}_sorted_reads.bam
	samtools index barcode${bc}_sorted_reads.bam
	mv barcode${bc}_sorted_reads.bam ../sorted_bam_files/
	mv barcode${bc}_sorted_reads.bam.bai ../sorted_bam_files/
done
#Output will be X sorted and indexed bam files, where X is the number of barcodes listed in barcode_number. There will be X .bam files and X .bam.bai files. These files will be saved to the new directory barcoded_location/sorted_bam_files. 
#Output of this step can be visualized in Tablet, or any other alignment viewing software. 

#STEP 5: Call SNPs for each barcoded individual against the reference genome. 
cd ${barcoded_location}/sorted_bam_files/
#The following line indexes the reference genome. We're using a different program here, so you have to do it again, even though you already did something similar when using BWA.
mv ${barcoded_location}/fastq_summary/${ref_genome} ${barcoded_location}/sorted_bam_files/
samtools faidx ${ref_genome}
#The following lines call SNPs for all your samples against the reference genome.
#The -d option allows you to specify the maximum per-sample depth. It's set to 100x right now, but if you suspect you might have deeper coverage than that, you can set it to as high as you'd like. We've done up to 1000000x.
#The -Ov option just specifies that your output should be a vcf file.
bcftools mpileup -b ${all_samples} -d 100 -f ${ref_genome} -o all_samples_mpileup.vcf -Ov 
#The -m option below is used to invoke the multiallelic caller.
#You can also optionally invoke the -V option below, which allows you to select a type of variant to skip (either indels or snps). We often chose to ignore indels.
bcftools call -Ov -m -o all_samples_called.vcf -v all_samples_mpileup.vcf
#Output will be two vcf files, one named all_samples_mpileup.vcf and one named all_samples_called.vcf. The one you want to work with is all_samples_called.vcf.

#STEP 6: Estimate Fst using vcftools
#The following code only includes provisions for two populations, but you could easily add more if you have individuals belonging to more than two populations. Just add another --weir-fst-pop option, and in the input section above, add pop_3 as a variable with its respective path.
vcftools --vcf all_samples_called.vcf --weir-fst-pop ${pop_1} --weir-fst-pop ${pop_2} --out pop_1_to_pop_2_fst.txt 
#Output will be a text file containing Fst values for sites that vary against the reference genome.

#Citations
#Danecek P, Schiffels S, Durbin R. 2014. Multiallelic calling model in bcftools (-m).
#Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, Handsaker RE, Lunter G, Marth GT, Sherry ST, McVean G, Durbin R, 1000 Genomes Project Analysis Group. 2011. The variant call format and VCFtools. Bioinformatics. 27(15): 2156-2158.
#Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, 1000 Genome Project Data Processing Subgroup. 2009. The sequence alignment/map (SAM) format and SAMtools. Bioinformatics. 25(16): 2078-2079.
#Li H. 2010. Mathematical notes on SAMtools algorithms.
#Li H. 2011. A statistical framework for SNP calling, mutation discovery, association mapping, and population genetical parameter estimation from sequencing data. Bioinformatics. 27(21): 2987-2993.
#Li H. 2013. Aligning sequence reads, clone sequences, and assembly contigs with BWA-MEM. arXiv:1303.3997