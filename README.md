# miRNA workflow

This workflow is an adaptation to run the miRNA pipeline developed by Pratibha Potla for Dr. Amanda Ali (MAR-13-2020):

## How to run this workflow

### 1) Setting up the environment

A) If you're running the pipeline in your computer or server, install or load the dependecies:
	bcl2fastq, gcc/9.3.0, samtools/1.13, bowtie/1.3.0, bedtools/2.29.2, python3.8, cutadapt-3.4

b) If you're running the pipeline in Compute Canada clusters create a miRNA_workflow environment:
      
	Module load python3.8
	module load scipy-stack
	virtualenv --no-download miRNA_workflow
      
	# Activate virtual environment
	source cutadapt/bin/activate
	# upgrade pip
	pip install --no-index --upgrade pip
	
	# install cutadapt (cutadapt/3.4)
	pip install cutadapt --no-index

	# deactivate environment
	deactivate


## Running final-pipeline_new.sh

### A) PREPARATION

#### If you're running the pipeline in your computer or server, follow instructions by Potla (2020):

	1) Download  both mature and hairpin fasta file from http://www.mirbase.org/ftp.shtml and extract only human (hsa) sequences from each one of them

	2) Create bowtie1 index files for both these fasta files, where files with extension .ebwt will be created

	3) Download latest human reference genome version GRCh38 from Ensembl or NCBI in FASTA format having chr1, chr2, etc all in one file

	4) Create bowtie1 index files for reference.fasta file, where files with extension .ebwt will be created

	5) If sample filename is 'Bone1B_S2_R1_001.fastq.gz', then rename it to 'Bone1B.fastq.gz' for convenience in using in scripts

	6) Download and install latest versions of these softwares: bcl2fastq, fastqc, umi_tools, cutadapt, bowtie1, samtools, and bedtools (tagBAM script only)

	7) Keep the file 'hsa-genome-miRBase22v-onlymiRNAs-convforTagBAM.bed' under the same directory where miRBase index files are present



#### If you're running the pipeline in Compute Canada clusters:

	1) Required databases are ready in the following paths:
		Human sequences (hsa) of both mature and hairpin fasta file (from http://www.mirbase.org/ftp.shtml): /lustre03/project/6004655/COMMUN/data/miRNA_databases
		
		Bowtie indexes: /lustre03/project/6004655/COMMUN/data/miRNA_databases

		Annotation File path: "/lustre03/project/6004655/COMMUN/data/miRNA_databases/hsa.gff3"
	
		bowtie1 index of the human reference genome version GRCh38 from Ensembl: /cvmfs/ref.mugqic/genomes/species/Homo_sapiens.GRCh38/genome/bowtie_index/
	
	2) You must modify the line "activate miRNA_workflow environment" in the final-pipeline_new_noUMIs.sh script according to the path you used.
	The default is to activate the preset miRNA_worflow environment in this location:
	/home/p1106294/miRNA_workflow

	3) Pipeline script "final-pipeline_new_noUMIs.sh" activates the environment and loads fastqc, bowtie1, samtools and bedtools.
	
	4) The file: 'hsa-genome-miRBase22v-onlymiRNAs-convforTagBAM.bed' must be under the same directory where miRBase index files are present.
	
	5) All bam files must be in the same directory. If there are fastq files, all fastq files must be in the same directory and, step
  "Convert Bam file to fastq" script chunk must be commented (add "#" and the beginning of each line in the the script "final-pipeline_new.sh").
  
  	6) If sample filename is 'Bone1B_S2_R1_001.fastq.gz', then rename it to 'Bone1B.fastq.gz' for convenience in using in scripts

	7) Create a new directory called "analysis" to store results


### B) EXECUTION:

1. Run the script "final-pipeline_new_noUMIs.sh" with this arguments :

	1. Input directory
	2. Output directory
	3. Readsets txt file
	4. Library type: "paired"  or "single"
	5. Sample map: when there are several read groups per sample, the sample map is a txt file mapping samples to readsets
	6. Path to human Mirbase data base index
	7. Path to human reference genome
	8. Path to human Mir Annotated gff3 file
	9. gencode.v39.chr_patch_hapl_scaff.annotation.gtf.gz file
	10. 3`adapter
	11. 5`adapter


	Example:

	bash final-pipeline_new_noUMIs.sh ~/miRNA_project/data ~/miRNA_project/analysis sample_list.txt "paired" RGsampleMap.txt ~/miRNA_project/miRNA_databases/mature_index \
	/cvmfs/ref.mugqic/genomes/species/Homo_sapiens.GRCh38/genome/bowtie_index/ ~/miRNA_project/miRNA_databases/hsa.gff3 "TGGAATTCTCGGGTGCCAAGG"

2. When 'final-pipeline_new_noUMIs.sh' has finished running, copy scripts 'get-count_one-sample.sh' and "get-count_job.sh" to analysis directory. 
   Copy a "sample_list.txt" file to analysis directory too.
   "get-count_job.sh" script launches an array job to count alinged to reference genome miRNAs in all samples in parallel. Launch this script with:

    sbatch get-count_job.sh

3. Run python custom script "merge_miRNA_counts_files.py". Arguments:

	1. Sample list
	2. Path to mirbase mature miRNA counts
	3. Path to aligned to reference genome miRNA counts
	4. Output directory

	Example:

	module load python/3.9
	module load scipy-stack
	python merge_miRNA_counts_files.py "sample_list.txt" "../analysis/maturemiRNAcounts" "../analysis/genome-miRNAcounts" "../analysis"





   

   





