# RNA-Seq pipeline for yeast.

#This is the bash script used to complete the RNA-Seq analysis from quality check to quantification.
#The script was compiled and edited in Notepad++

#On running this script for a set of yeast RNA-seq data against the reference yeast genome, we get the quantified gene expression data as the final output. 

#!/bin/bash

# We initiate our RNA-Seq pipeline
echo "Start RNA-Seq Pipeline"

# Set the working directory
cd /hara

#Step 1: Get the files from sources using wget and place in the working directory 

#In our analysis, we use the paired RNA-seq files which are in fastq format as submitted by to EBI (European Nucleotide Archive) 
#We can download the files directly from the web browser or from within the command prompt using the wget command.

#File 1 (ERR3772427_1.fastq.gz) is the forward reads of the sample dataset.
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR377/007/ERR3772427/ERR3772427_1.fastq.gz

#Unzip the file into the working directory
gzip -d Project/ERR3772427_1.fastq

#File 2 (ERR3772427_2.fastq.gz) is the reverse reads of the sample dataset 
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR377/007/ERR3772427/ERR3772427_2.fastq.gz

#Unzip the file into the working directory
gzip -d Project/ERR3772427_2.fastq

#Now, our .fastq files are ready to have their quality of reads checked with the first tool which is FASTQC

#Step 2: Run fastqc on the RNA-Seq foward and reverse files to check for quality of reads.

echo "Starting FASTQC for RNA-Seq dataset"
echo "------------------------------------------------------------------------------------------------"

fastqc Project/ena_files/ERR3772427_1.fq Project/ena_files/ERR3772427_2.fq -o Project/ena_files

echo "------------------------------------------------------------------------------------------------"
echo "FASTQC finished running" 

#Once FASTQ finishes running, we get the Quality Check report in the form of .html files which can be opened in a web browser. 
#The files are attached in the Project folder which has been submitted as a part of the project. 

#They are named ERR3772427_1.fastq and ERR3772427_2.fastq respectively. 

#On analysing the report, we can see that trimming of the poor reads must be carried out, as some reads are in the yellow zone.
#Hence, we input the forward and reverse files into Trimmomatic and get output 4 files 
#(Unpaired Forward, Paired Forward, Unpaired Reverse and Paired Reverse). 

echo "------------------------------------------------------------------------------------------------"

#Step 3: Our second tool is Trimmomatic which is used to trim the poor reads and obtain 4 output files 

echo "------------------------------------------------------------------------------------------------"
echo "Trimmomatic starting for RNA-Seq dataset"
echo "------------------------------------------------------------------------------------------------"

#We run trimmomatic for the forward and reverse files with the set conditions.
#We use illumina adapters pre-configured in trimmomatic for getting the paired outputs. SInce our input data is sequenced using Illumina, we use the same adapters.
#We set the minimum read length to 25 which can be changed to suit our data if necessary. 

java -jar ~/Trimmomatic-0.39/Trimmomatic-0.39/trimmomatic-0.39.jar PE Project/ena_files/ERR3772427/ERR3772427_1.fastq Project/ena_files/ERR3772427/ERR3772427_2.fastq ERR3772427_1_forward_paired.fq.gz ERR3772427_1_forward_unpaired.fq.gz ERR3772427_2_reverse_paired.fq.gz ERR3772427_2_reverse_unpaired.fq.gz ILLUMINACLIP:Trimmomatic-0.39/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25

echo "------------------------------------------------------------------------------------------------"
echo "Trimmomatic finished running!"

# Output file names called ERR3772427_1_forward_paired.fq.gz , ERR3772427_1_forward_unpaired.fq.gz ,ERR3772427_2_reverse_paired.fq.gz and ERR3772427_2_reverse_unpaired.fq.gz

# We will use ERR3772427_1_forward_paired.fq.gz and ERR3772427_2_reverse_paired.fq.gz for our further analyses, as the reads of interest are only the paired reads.


#Step 4: Run fastqc again for the PAIRED files to check the results of the trimming.
 
#Now, we run FASTQC again to check the quality of the new reads.

echo "Starting FASTQC for trimmed and paired files"
echo "------------------------------------------------------------------------------------------------"

fastqc Project/ena_files/ERR3772427_1_forward_paired.fq Project/ena_files/ERR3772427_2_reverse_paired.fq -o Project/ena_files

echo "------------------------------------------------------------------------------------------------"
echo "FASTC for trimmed and paired files finished running"

#Once FASTQ finishes running, we get the new Quality Check report .html files which can be opened in the web browser. 

#The files are attached in the Project folder which has been submitted as a part of the assignment. 
#They are named ERR3772427_1_forward_paired_fastqc.html and ERR3772427_2_reverse_paired_fastqc.html respectively. 

#On analysing the report, we can see that the reads are of significantly better quality, indicated that the trimming using Trimmomatic was successful. 

#For the next step, we use our third tool in the pipeline which is HISAT2. 
#We use this tool for the alignment of the paired sequences to the reference yeast genome.

#Step 5: Alignment of the reads using HISAT2 to the reference yeast genome

echo "------------------------------------------------------------------------------------------------"
echo "Starting HISAT2 for alignment of reads to reference genome"
echo "------------------------------------------------------------------------------------------------"

# Our objective is to input the paired files and align it against the yeast genome to get an output SAM file called yeastaligned.sam

#Before running our tool, we must download the yeast reference genome in INDEX form. 
#It can be found on the HISAT2 website and retreived using wget or downloaded manually from the website.

wget  https://cloud.biohpc.swmed.edu/index.php/s/JRSoKHD5cHfpCFE/download

# The download file must then be unzipped and unpacked. 

gzip -d r64.tar.gz 
tar -xf r64.tar 

#Now that the index reference genome file is ready, we can run HISAT2 and align the paired sequence files.

hisat2 -q -x hisat2-2.2.1/r64/genome -1 Project/ena_files/ERR3772427_1_forward_paired.fq -2 Project/ena_files/ERR3772427_2_reverse_paired.fq  -S hisat2-2.2.1/yeastaligned.sam

echo "------------------------------------------------------------------------------------------------"
echo "HISAT2 successfully run, sequences aligned"
echo "------------------------------------------------------------------------------------------------"

#Once we obtain our .SAM file, we are required to convert it into a sorted .BAM file to facilitate the quantification analysis.

#The conversion is done using the following command, which will give an unsorted output .BAM file yeastaligned.bam 

samtools view -S -b hisat2-2.2.1/yeastaligned.sam > hisat2-2.2.1/yeastaligned.bam

#However, the .BAM file must still be converted to a sorted .BAM file, and that is done by giving the following command.

samtools sort hisat2-2.2.1/yeastaligned.bam -o hisat2-2.2.1/yeastaligned_sorted.bam

#Now, our sorted .BAM file is ready for quantification. 

#Step 6: We can finally use our fourth tool of the pipeline featureCounts, which will quantify the gene expression. 

#However, before quantification we require a .gtf file which contains the coordinates of genomic features of the data. 
#It is required for the accurate and precise quantification of reads. 

# We can get the yeast gtf file from the Ensembl server using wget or downloading directly through the browser.

wget  http://ftp.ensembl.org/pub/release-107/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.107.gtf.gz

# The download file must then be unzipped 
gzip -d Saccharomyces_cerevisiae.R64-1-1.107.gtf.gz

#Now, our gtf file is ready to facilitate the quantification. 

echo "------------------------------------------------------------------------------------------------"
echo "Starting quantification of gene expression using featureCounts"
echo "------------------------------------------------------------------------------------------------"

featureCounts -S 2 -a Saccharomyces_cerevisiae.R64-1-1.107.gtf -o Project/yeast_featurecounts.txt hisat2-2.2.1/yeastaligned_sorted.bam

echo "Quantification complete, summary text files created"
echo "------------------------------------------------------------------------------------------------"

#The quantification is now complete and the data may be accessed in the yeast_featurecounts.txt file.
#In addition, we also get a summary file which represents the highlights of the dataset including factors such as ..


#We can now determine the top 5 highest expressed genes in our given yeast dataset. 
#We first slice the columns of interest, which are the GENEID and their counts, and then sort them in descending order. We then show the top 5 entries in the file.
echo "The 5 highest expressed genes are:\n"

cat Project/yeast_featurecounts.txt |cut -f 1,7 | sort -k 2 -r -g | head -5

echo "------------------------------------------------------------------------------------------------"

#The output of the above command is 
#YER065C 227627
#YKR097W 169503
#YAL054C 126319
#YOR374W 123019
#YER024W 115222

#We can also find the number of genes with no expression (Zero count). Again, we extract columns of interest GENEID and their counts, and then extract the lines with 0. 
#Then count the number of lines with zero. 

echo "The number of genes with no expression is:\n"

cat Project/yeast_featurecounts.txt | cut -f 1,7 | tail -n +2 | grep -w 0 | wc -l

echo "------------------------------------------------------------------------------------------------"

#The output of the above command is 
# 497





