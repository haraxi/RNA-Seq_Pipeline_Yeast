# RNA-Seq Pipeline Yeast

## **Course Name** 
Practical Computing for Scientists (03-701)
Carnegie Mellon University, Pittsburgh, PA

## **Overview**
This is my implementation of an RNA-seq pipeline using Yeast RNA-seq data from the European Nucleotide Archive (Project Accession: [PRJEB35903](https://www.ebi.ac.uk/ena/browser/view/PRJEB35903)). 
Here, I wrote a Bash script to carry out integral steps in RNA-seq analysis including quality check, trimming, alignment and quantification, effictively creating a pipeline to carry out these tasks.
Each step utilised a different linux-based bioinformatics tool to carry out the specific task in the pipeline. 

1. Quality Check: [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. Trimming: [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
3. Alignment to reference genome: [HISAT2](http://daehwankimlab.github.io/hisat2/)
4. Quantification: [FeatureCounts](https://subread.sourceforge.net/featureCounts.html)

## **Results**
On running the complete script, I was able to obtain an overall alignment of 95.91% to the reference yeast genome, and quantify the gene expression levels. On quantification, I identified the number of genes with no expression in the control sample as well as the genes with the highest expression levels. A detailed explanation of the project workflow and results can be found in [the project report](FinalProjectReport.pdf).

Future work on this project will include running the pipeline on a treatment sample as well, and carrying out a comparative analysis of gene expression in the control and sample using differential expression analysis.
