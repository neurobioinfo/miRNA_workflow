#!/bin/sh

# exit when any command fails
# set -e

# keep track of the last executed command
# trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
# trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

args_number=$#
echo $args_number
if [[ $args_number -eq 0 ]]; then
    echo "Usage: "
    echo "bash final-pipeline_new.sh <input dir> \\"
    echo "< output dir > \\"
    echo "< readset list > \\"
    echo "< single or paired > \\"
    echo "< sample map file > \\"
    echo "< path to human Mirbase Index > \\"
    echo "< path to human Reference genome Index > \\"
    echo "< path to human Mir Annotated gg3 file > \\"
    echo "< human Mir TagBam Bed file > \\"
    echo "< adapter 3' adapter 5' >"
    echo
    exit 1
fi

echo "Activating miRNA_workflow environment"
echo "Uncomment if you're running the pipeline in Compute Canada Clusters"
echo "Comment if you're using you're own installation"
module load scipy-stack
source ~/miRNA_workflow/bin/activate

echo "Loading modules"
module load fastqc/0.11.9
module load StdEnv/2020 gcc/9.3.0 samtools/1.13
module load bowtie/1.3.0
module load bedtools/2.29.2

# Passing arguments
inputdir=$1
outputdir=$2
readSetlist=$3
libraryLayout=$4
sampleMap=$5
humanMirbaseIndex=$6
humanRefIndex=$7
humanMirAnnotated=$8
humanMirTagBamBed=$9
adapter3=${10}
adapter5=${11}



mkdir -p $inputdir/QCReports
mkdir -p $outputdir/maturemiRNAcounts
mkdir -p $outputdir/genome-miRNAcounts
mkdir -p $outputdir/tmp
export TMPDIR=$outputdir/tmp

# For each read set:
while read i;
do
    #echo "Convert Bam file to fastq"
    #samtools collate -u -O $inputdir/$i.bam | \
    #samtools fastq -@ 16 -1 $inputdir/${i}_R1.fastq -2 $inputdir/${i}_R2.fastq -0 /dev/null -s /dev/null -n
    #echo "Bam file to fastq done"
    #echo
    
    echo "Running FastQC on raw sample $i now"
    if [[ $libraryLayout == "paired" ]]; then
        fastqc -t 16 -o $inputdir/QCReports $inputdir/$i"_R1.fastq" $inputdir/$i"_R2.fastq"
    else
        fastqc -t 16 -o $inputdir/QCReports $inputdir/$i".fastq"
    fi
    mkdir $outputdir/$i
    echo

    # UMI extraction
    #echo "Sample $i: Running of UMI extraction with regex on raw reads"
    echo "UMI extraction disabled for this miRNA workflow version"
    echo
    #umi_tools extract --stdin=$inputdir/$i".fastq" --log=$outputdir/$i/$i"-UMIextraction-fromrawreads.log" 
    #   --stdout=$outputdir/$i/$i"-directUMIextracted.fastq" --extract-method=regex 
    #   --bc-pattern='.+(?P<discard_1>AACTGTAGGCACCATCAAT){s<=2}(?P<umi_1>.{12})(?P<discard_2>.+)'

    # Remove too short and too long reads with cutadapt
    echo "Sample $i: Applying a read length filter to remove too short and long reads"
    if [[ $libraryLayout == "paired" ]]; then
        cutadapt -j 16 -q 10 -a ${adapter3} -A ${adapter5} --minimum-length=15 --maximum-length=32 \
        -o $outputdir/$i/$i"_R1-min18max30L.fastq" \
        -p $outputdir/$i/$i"_R2-min18max30L.fastq" \
        $inputdir/$i"_R1.fastq" $inputdir/$i"_R2.fastq" 2> $outputdir/$i/$i"-readlengthfilter-cutadapt.log"
        echo "adapter removing and trimming done"
    else
        echo $libraryLayout $adapter3
        cutadapt -j 16 -q 10 -a ${adapter3} --minimum-length=15 --maximum-length=32 \
        -o ${outputdir}/${i}/${i}"_min18max30L.fastq" ${inputdir}/${i}.fastq 2> ${outputdir}/${i}/${i}-readlengthfilter-cutadapt.log
        echo "adapter removing and trimming done"
    fi
    echo

    # Align to MirBase database with Bowtie
    echo "Running bowtie to mirbase mature seq"
    # Paired end option
    if [[ $libraryLayout == "paired" ]]; then
        bowtie -q -n 1 -l 20 --allow-contain -k 1 --threads 16  -S -x $humanMirbaseIndex/"hsa_mature" \
        -1 $outputdir/$i/${i}"_R1-min18max30L.fastq" -2 $outputdir/$i/${i}"_R2-min18max30L.fastq" \
        --un $outputdir/$i/$i"-maturemiRNA-unalignedReads-bowtie1-beststratam1.fastq.gz" \
        $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-beststratam1.sam" 2> $outputdir/$i/$i"-bowtie-maturemiRNA-beststratam1.log"
        echo "Bowtie alignment done"

        echo "Sorting the BAM file"
        echo "Sorting by name to fix mate, When there are not UMIs, no other steps are needed"
        samtools sort -n -@ 16 $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-beststratam1.sam" > $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-beststratam1.bam"
        samtools fixmate -m $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-beststratam1.bam" $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-fixmate.bam"
        samtools sort -n -@ 16 $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-fixmate.bam" > $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-nsorted.bam"
        rm -rf $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-beststratam1.sam"
        rm -rf $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-fixmate.bam"
        echo "Sorting and indexing done"
        echo

        echo "Running of bowtie from unaligned mature miRNA reads to genome as reference (Exiqon style)"
        bowtie -q -n 1 -l 20 --allow-contain -k 1 --threads 16  -S -x $humanRefIndex/"Homo_sapiens.GRCh38" \
        -1 $outputdir/$i/${i}"-maturemiRNA-unalignedReads-bowtie1-beststratam1.fastq_1.gz" \
        -2 $outputdir/$i/${i}"-maturemiRNA-unalignedReads-bowtie1-beststratam1.fastq_2.gz" \
        --al $outputdir/$i/${i}"-genomeaftermiRNA-alignedReads-bowtie1.fastq" -S $outputdir/$i/$i"-genomeaftermiRNA-alignedReads-bowtie1.sam" 2> $outputdir/$i/$i"-bowtie-genomeaftermiRNA.log"
        echo "Bowtie of unaligned mature miRNA reads to genome reference (Exiqon style) done"
        echo

        echo "Sorting and indexing for downstream analyses"
        samtools sort -n -@ 16 $outputdir/$i/$i"-genomeaftermiRNA-alignedReads-bowtie1.sam" > $outputdir/$i/$i"-genomeaftermiRNA-aligned-bowtie1.bam"
        samtools fixmate -m $outputdir/$i/$i"-genomeaftermiRNA-aligned-bowtie1.bam" $outputdir/$i/$i"-genomeaftermiRNA-aligned-bowtie1-fixmate.bam"
        samtools sort -n -@ 16 $outputdir/$i/$i"-genomeaftermiRNA-aligned-bowtie1-fixmate.bam" > $outputdir/$i/$i"-genomeaftermiRNA-aligned-bowtie1-nsorted.bam"
        rm -rf $outputdir/$i/$i"-genomeaftermiRNA-alignedReads-bowtie1.sam"
        rm -rf $outputdir/$i/$i"-genomeaftermiRNA-aligned-bowtie1-fixmate.bam"
        echo

    # Single end option
    else
        # Alignment to MirBase database
        bowtie -q -v 1 -l 30 -a --best --strata --threads 16  -S -x $humanMirbaseIndex/"hsa_mature" \
        $outputdir/$i/${i}"_min18max30L.fastq" \
        --un $outputdir/$i/$i"-maturemiRNA-unalignedReads-bowtie1-beststratam1.fastq.gz" $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-beststratam1.sam" \
        2>$outputdir/$i/$i"-bowtie-maturemiRNA-beststratam1.log"
        echo "Bowtie alignment done"
        echo

        echo "Sorting by name. When there are not UMIs, no other steps are needed"
        samtools sort -n -@ 16 $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-beststratam1.sam" > $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-nsorted.bam"
        echo "Sorting done"
        echo

        #Alignment to reference genome h38
        echo "Running of bowtie from unaligned mature miRNA reads to genome as reference (Exiqon style)"
        bowtie -q -v 1 -l 30 -a --best --strata --threads 16  -S -x $humanRefIndex/"Homo_sapiens.GRCh38" \
        $outputdir/$i/$i"-maturemiRNA-unalignedReads-bowtie1-beststratam1.fastq.gz" \
        --al $outputdir/$i/${i}"-genomeaftermiRNA-alignedReads-bowtie1.fastq" -S $outputdir/$i/$i"-genomeaftermiRNA-alignedReads-bowtie1.sam" 2> $outputdir/$i/$i"-bowtie-genomeaftermiRNA.log"
        echo "Bowtie of unaligned mature miRNA reads to genome reference (Exiqon style) done"
        echo 

        # Sorting by  name
        echo "Sorting for downstream analyses"
        samtools sort -n -@ 16 $outputdir/$i/$i"-genomeaftermiRNA-alignedReads-bowtie1.sam" > $outputdir/$i/$i"-genomeaftermiRNA-aligned-bowtie1-nsorted.bam"
        echo "samtools sorting by name done"

    fi
    #break
    done < $readSetlist


#For multiple readsets per sample stop here to merge. then continue
# By the moment, working only for single end
echo "Merging readsets"
while read line; do
    #read Sample-readset map
    nline=$(echo $line|wc -w)
    echo $nline
    nRG=$((nline - 1))
    sample=$(echo $line|awk '{print $1}')
    echo "Merge readsets of "$sample
    mkdir $outputdir/$sample

    # If there is only one readset
    if [[ $nRG == 0 ]]; then
        echo "readset number "$((nRG + 1))
        mv $outputdir/$sample/$sample"-maturemiRNA-aligned-bowtie1-nsorted.bam" $outputdir/$sample/$sample"-maturemiRNA-aligned-bowtie1-merged.bam"
        mv $outputdir/$sample/$sample"-genomeaftermiRNA-aligned-bowtie1-nsorted.bam" $outputdir/$sample/$sample"-genomeaftermiRNA-aligned-bowtie1-merged.bam"

    else
    #For multiple readsets
        touch bamfiles_mature
        touch bamfiles_genome
        RGs=$(echo $line|awk '{for (i=2; i<=NF; i++) print $i}')
        for RG in $RGs;do
            echo $outputdir/$RG/$RG"-maturemiRNA-aligned-bowtie1-nsorted.bam" >> bamfiles_mature
            echo $outputdir/$RG/$RG"-genomeaftermiRNA-aligned-bowtie1-nsorted.bam" >> bamfiles_genome
            mv $outputdir/$RG/*.log $outputdir/$sample/
            mv $outputdir/$RG/${RG}"-genomeaftermiRNA-alignedReads-bowtie1.fastq" $outputdir/$sample/
            done 
        echo "Merging mature alingment files"
        samtools merge -n -f --no-PG -@ 16 -b bamfiles_mature -o $outputdir/$sample/$sample"-maturemiRNA-aligned-bowtie1-merged.bam"
        echo "Merging reference genome hg38 alingment files"
        samtools merge -n -f --no-PG -@ 16 -b bamfiles_genome -o $outputdir/$sample/$sample"-genomeaftermiRNA-aligned-bowtie1-merged.bam"
        for RG in $RGs;do
            rm -r $outputdir/$RG
            done
        rm bamfiles_mature
        rm bamfiles_genome
    fi
    echo

    # Sorting by coordinates and indexing
    echo "Sorting by coordinates and indexing maturemiRNA-aligned reads"
    # MirBase alignments
    samtools sort -@ 16 $outputdir/$sample/$sample"-maturemiRNA-aligned-bowtie1-merged.bam" > $outputdir/$sample/$sample"-maturemiRNA-aligned-bowtie1-CoordSort.bam"
    samtools index -@ 16 $outputdir/$sample/$sample"-maturemiRNA-aligned-bowtie1-CoordSort.bam"
    rm -rf $outputdir/$sample/$sample"-maturemiRNA-aligned-bowtie1-merged.bam"
    #rm -rf $outputdir/$sample/$sample"-maturemiRNA-aligned-bowtie1-fixmate.bam"
    echo

    # Reference genome alignments
    #echo "Sorting and indexing genome-aligned readsfor downstream analyses"
    samtools sort -@ 16 $outputdir/$sample/$sample"-genomeaftermiRNA-aligned-bowtie1-merged.bam" > $outputdir/$sample/$sample"-genomeaftermiRNA-aligned-bowtie1-CoordSort.bam"
    samtools index -@ 16 $outputdir/$sample/$sample"-genomeaftermiRNA-aligned-bowtie1-CoordSort.bam"
    rm -rf $outputdir/$sample/$sample"-genomeaftermiRNA-aligned-bowtie1-merged.bam" 
    #rm -rf $outputdir/$i/$i"-genomeaftermiRNA-aligned-bowtie1-fixmate.bam"
    echo "Sorting and indexing done"
    echo

    echo "Counting step"
    #echo "Couting step without UMI's will not use Umi tools. Continue counting with samtools idxstats"
    #umi_tools count --method=unique --per-contig -I $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-beststratam1.bam" \
    #-L $outputdir/$i/$i"_counts-uniquemethod-maturemiRNA.log" -S $outputdir/$i/$i"_counts-finaloutput-uniquemethod-maturemiRNA.txt"
    #echo

    #echo "Deduplicating the aligned BAM"
    #echo "Removal of duplicated read isn't recomended for RNAseq data without UMIs"
    #umi_tools dedup --method=unique -I $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-beststratam1.bam" \
    #-S $outputdir/$i/$i"_deduplicated-matureMirna-uniquemethod-beststratam1.bam" -L $outputdir/$i/$i"-deduplicate-matureMirna-uniquemethod-beststratam1.log"
    #echo
    #echo "Indexing the BAM output for finding counts"
    #echo "Not needed for data without UMIS"
    #samtools index $outputdir/$i/$i"_deduplicated-matureMirna-uniquemethod-beststratam1.bam"
    #echo
    echo

    echo "MirBase miRNA counting for sample "${sample}
    samtools idxstats -@ 16 $outputdir/$sample/$sample"-maturemiRNA-aligned-bowtie1-CoordSort.bam" | cut -f1,3 - | sed "1s/^/miRNA\t${i}-miRNAcount\n/" - > $outputdir/maturemiRNAcounts/${sample}"-maturemiRNA-counts.txt"
    echo "Done"
    echo


    ## This isn't necessary for data with no UMIs
    #echo "Deduplicating the BAM file before counting through custom scripts"
    ##umi_tools dedup --method=unique -I $outputdir/$i/$i"-genomeaftermiRNA-aligned-bowtie1-beststratam1.bam" -S $outputdir/$i/$i"_deduplicated-genomeaftermiRNA-uniquemethod.bam" -L $outputdir/$i/$i"-deduplicate-genomeaftermiRNA-uniquemethod.log"
    #samtools markdup -l 30 -r -T /home/p1106294/scratch/2021-05-25.sara/analysis/tmp/${i} -f $outputdir/$i/${i}".duplicate-align-genome.stats" -m s \
    #--include-fails -@ 16 $outputdir/$i/$i"-genomeaftermiRNA-aligned-bowtie1.bam" $outputdir/$i/$i"-maturemiRNA-genomealigned-bowtie1-markdup.bam"
    #echo
    #echo "Index the deduplicated BAM file"
    #samtools index $outputdir/$i/$i"-maturemiRNA-genomealigned-bowtie1-markdup.bam"
    #echo

    echo "Genome aligned miRNA counting for sample "${sample}
    echo "Actual filtering sam step only for overlapping microRNA locations"
    tagBam -i $outputdir/$sample/$sample"-genomeaftermiRNA-aligned-bowtie1-CoordSort.bam" -files $humanMirTagBamBed \
    -names -tag XQ > $outputdir/$sample/$sample"-tagged.bam"
    echo

    #break
    done < $sampleMap

echo "Completed analysis for whole run!!"
exit 0
