#!/bin/sh

# exit when any command fails
# set -e

# keep track of the last executed command
# trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
# trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT



echo "Activating miRNA_workflow environment"
module load python/3.8.10
module load scipy-stack
source ~/miRNA_workflow/bin/activate

echo "Loading modules"
module load fastqc/0.11.9
module load StdEnv/2020 gcc/9.3.0 samtools/1.13
module load bowtie/1.3.0
module load bedtools/2.29.2

# Passing arguments
inputdir="$1"
echo inputdir $inputdir
outputdir="$2"
echo outputdir $outputdir
libraryLayout="$3"
echo libraryLayout $libraryLayout
humanMirbaseIndex="$4"
echo humanMirbaseIndex $humanMirbaseIndex
humanRefIndex="$5"
echo humanRefIndex $humanRefIndex
humanMirAnnotated="$6"
echo humanMirAnnotated $humanMirAnnotated
adapter3="$7"
echo adapter3 $adapter3
adapter5="$8"
echo adapter5 $adapter5
sampleRG=($9)
echo sampleRG
echo sampleRG "${sampleRG[@]}"

echo


humanMirTagBamBed=$humanMirbaseIndex"/hsa-genome-miRBase22v-onlymiRNAs-convforTagBAM.bed"

mkdir -p $inputdir/QCReports
mkdir -p $outputdir/maturemiRNAcounts
mkdir -p $outputdir/genome-miRNAcounts
mkdir -p $outputdir/tmp
export TMPDIR=$outputdir/tmp

# Read Sample-readset map
echo "spliting sample and RGs"
echo sampleRG "${sampleRG[@]}"
sample="${sampleRG[0]}"
echo sample $sample
unset sampleRG[0]

nRG=${#sampleRG[@]}
echo Number of readsets $nRG

# For read set:
for i in ${sampleRG[@]}; do

    echo Processing Read Set $i
    echo

    #echo "Convert Bam file to fastq"
    #samtools collate -u -O $inputdir/$i.bam | \
    #samtools fastq -@ 16 -1 $inputdir/${i}_R1.fastq -2 $inputdir/${i}_R2.fastq -0 /dev/null -s /dev/null -n
    #echo "Bam file to fastq done"
    #echo

    #echo "Running FastQC on raw sample $i now"
    #if [[ $libraryLayout == "paired" ]]; then
        #fastqc -t 16 -o $inputdir/QCReports $inputdir/$i"_R1.fastq" $inputdir/$i"_R2.fastq"
    #else
        #fastqc -t 16 -o $inputdir/QCReports $inputdir/$i".fastq"
    #fi

    mkdir $outputdir/$i
    echo

    # UMI extraction
    echo "Read set $i: Running of UMI extraction with regex on raw reads"
    #echo "Skiping UMI extraction for data with no UMIs"
    umi_tools extract --stdin=$inputdir/$i".fastq.gz" --log=$outputdir/$i/$i"-UMIextraction-fromrawreads.log" \
        --stdout=$outputdir/$i/$i"-directUMIextracted.fastq.gz" --extract-method=regex \
        --bc-pattern='.+(?P<discard_1>AACTGTAGGCACCATCAAT){s<=2}(?P<umi_1>.{12})(?P<discard_2>.+)'
    echo

    # Remove too short and too long reads with cutadapt
    echo "Sample $i: Applying a read length filter to remove too short and long reads"
    ## Setup for paired read with UMIs TO DO
    if [[ $libraryLayout == "paired" ]]; then
        cutadapt -j 16 -q 10 --minimum-length=15 --maximum-length=32 \
        -o $outputdir/$i/$i"_R1-min18max30L.fastq" \
        -p $outputdir/$i/$i"_R2-min18max30L.fastq" \
        $inputdir/$i"_R1.fastq" $inputdir/$i"_R2.fastq" 2> $outputdir/$i/$i"-readlengthfilter-cutadapt.log"
        echo "adapter removing and trimming done"
    else
        echo Library layout $libraryLayout 
        cutadapt -j 16 -q 10 --minimum-length=15 --maximum-length=32 \
        -o ${outputdir}/${i}/${i}"_min18max30L.fastq.gz" $outputdir/$i/$i"-directUMIextracted.fastq.gz" 2>&1 | tee ${outputdir}/${i}/${i}-readlengthfilter-cutadapt.log
        echo "adapter removing and trimming done"
    fi
    echo

    # Align to MirBase database with Bowtie
    echo "Running bowtie to mirbase mature seq"
    # Paired end option SETUP TO DO
    if [[ $libraryLayout == "paired" ]]; then
        bowtie -q -n 1 -l 20 --allow-contain -k 1 --threads 16  -S -x $humanMirbaseIndex/"hsa_mature" \
        -1 $outputdir/$i/${i}"_R1-min18max30L.fastq" -2 $outputdir/$i/${i}"_R2-min18max30L.fastq" \
        --un $outputdir/$i/$i"-maturemiRNA-unalignedReads-bowtie1-beststratam1.fastq.gz" \
        $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-beststratam1.sam" 2> $outputdir/$i/$i"-bowtie-maturemiRNA-beststratam1.log"
        echo "Bowtie alignment done"
        echo "Sorting the BAM file"
        echo "Sorting by name to fix mate, When there are not UMIs, no other steps are needed"
        samtools sort -n -@ 16 -T $outputdir/"tmp" $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-beststratam1.sam" > $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-beststratam1.bam"
        samtools fixmate -m $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-beststratam1.bam" $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-fixmate.bam"
        samtools sort -n -@ 16  -T $outputdir/"tmp" $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-fixmate.bam" > $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-nsorted.bam"
        rm -rf $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-beststratam1.sam"
        rm -rf $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-fixmate.bam"
        echo "Sorting and indexing done"

        echo "Running of bowtie from unaligned mature miRNA reads to genome as reference (Exiqon style)"
        bowtie -q -n 1 -l 20 --allow-contain -k 1 --threads 16  -S -x $humanRefIndex/"Homo_sapiens.GRCh38" \
        -1 $outputdir/$i/${i}"-maturemiRNA-unalignedReads-bowtie1-beststratam1.fastq_1.gz" \
        -2 $outputdir/$i/${i}"-maturemiRNA-unalignedReads-bowtie1-beststratam1.fastq_2.gz" \
        --al $outputdir/$i/${i}"-genomeaftermiRNA-alignedReads-bowtie1.fastq" -S $outputdir/$i/$i"-genomeaftermiRNA-alignedReads-bowtie1.sam" 2> $outputdir/$i/$i"-bowtie-genomeaftermiRNA.log"
        echo "Bowtie of unaligned mature miRNA reads to genome reference (Exiqon style) done"
        echo
        echo "Sorting and indexing for downstream analyses"
        samtools sort -n -@ 16 -T $outputdir/"tmp" $outputdir/$i/$i"-genomeaftermiRNA-alignedReads-bowtie1.sam" > $outputdir/$i/$i"-genomeaftermiRNA-aligned-bowtie1.bam"
        samtools fixmate -m $outputdir/$i/$i"-genomeaftermiRNA-aligned-bowtie1.bam" $outputdir/$i/$i"-genomeaftermiRNA-aligned-bowtie1-fixmate.bam"
        samtools sort -n -@ 16 -T $outputdir/"tmp" $outputdir/$i/$i"-genomeaftermiRNA-aligned-bowtie1-fixmate.bam" > $outputdir/$i/$i"-genomeaftermiRNA-aligned-bowtie1-nsorted.bam"
        rm -rf $outputdir/$i/$i"-genomeaftermiRNA-alignedReads-bowtie1.sam"
        rm -rf $outputdir/$i/$i"-genomeaftermiRNA-aligned-bowtie1-fixmate.bam"
        echo

    # Single end option
    else
        # Alignment to MirBase database
        echo "Alignment to MirBase database Single end library layout"
        bowtie -q -v 1 -l 30 -a --best --strata --threads 16  -S -x $humanMirbaseIndex/"hsa_mature" \
        $outputdir/$i/${i}"_min18max30L.fastq.gz" \
        --un $outputdir/$i/$i"-maturemiRNA-unalignedReads-bowtie1-beststratam1.fastq.gz" $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-beststratam1.sam" \
        2>$outputdir/$i/$i"-bowtie-maturemiRNA-beststratam1.log"
        echo "Bowtie alignment done"
        echo

        echo "Sorting by name"
        samtools sort -n -@ 16 -T $outputdir/"tmp" $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-beststratam1.sam" > $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-nsorted.bam"
        echo "Sorting done"
        echo

        #Alignment to reference genome h38
        echo "Running of bowtie from unaligned mature miRNA reads to genome as reference (Exiqon style)"
        bowtie -q -v 1 -l 30 -a --best --strata --threads 16  -S -x $humanRefIndex/"Homo_sapiens.GRCh38" \
        $outputdir/$i/$i"-maturemiRNA-unalignedReads-bowtie1-beststratam1.fastq.gz" \
        --al $outputdir/$i/${i}"-genomeaftermiRNA-alignedReads-bowtie1.fastq.gz" \
        -S $outputdir/$i/$i"-genomeaftermiRNA-alignedReads-bowtie1.sam" 2> $outputdir/$i/$i"-bowtie-genomeaftermiRNA.log"
        echo "Bowtie of unaligned mature miRNA reads to genome reference (Exiqon style) done"
        echo 

        # Sorting by  name
        echo "Sorting for downstream analyses"
        samtools sort -n -@ 16 -T $outputdir/"tmp" $outputdir/$i/$i"-genomeaftermiRNA-alignedReads-bowtie1.sam" > $outputdir/$i/$i"-genomeaftermiRNA-aligned-bowtie1-nsorted.bam"
        echo "samtools sorting by name done"
    fi
    #break
done

#For multiple readsets per sample stop here to merge. then continue
# By the moment, working only for single end


#read Sample-readset map
echo "Merge readsets of sample "$sample
mkdir $outputdir/$sample
echo

# If there is only one readset
if [[ $nRG == 0 ]]; then
    echo "readset number "$((nRG + 1))
    mv $outputdir/$sample/$sample"-maturemiRNA-aligned-bowtie1-nsorted.bam" $outputdir/$sample/$sample"-maturemiRNA-aligned-bowtie1-merged.bam"
    mv $outputdir/$sample/$sample"-genomeaftermiRNA-aligned-bowtie1-nsorted.bam" $outputdir/$sample/$sample"-genomeaftermiRNA-aligned-bowtie1-merged.bam"

else
#For multiple readsets
    echo "Sample "$sample" has "$nRG" read sets "
    echo
    for RG in ${sampleRG[@]};do
        echo $outputdir/$RG/$RG"-maturemiRNA-aligned-bowtie1-nsorted.bam" >> ${sample}_bamfiles_mature
        echo $outputdir/$RG/$RG"-genomeaftermiRNA-aligned-bowtie1-nsorted.bam" >> ${sample}_bamfiles_genome
        mv $outputdir/$RG/*.log $outputdir/$sample/
        mv $outputdir/$RG/${RG}"-genomeaftermiRNA-alignedReads-bowtie1.fastq.gz" $outputdir/$sample/
    done

    echo "Merging maturemiRNA alingment files"
    samtools merge -n -f --no-PG -@ 16 -b ${sample}_bamfiles_mature -o $outputdir/$sample/$sample"-maturemiRNA-aligned-bowtie1-merged.bam"
    echo Done
    echo
    echo "Merging reference genome hg38 alingment files"
    samtools merge -n -f --no-PG -@ 16 -b ${sample}_bamfiles_genome -o $outputdir/$sample/$sample"-genomeaftermiRNA-aligned-bowtie1-merged.bam"
    echo Done
    echo
    #for RG in ${sampleRG[@]};do
    #    rm -r $outputdir/$RG
    #done
    rm ${sample}_bamfiles_mature
    rm ${sample}_bamfiles_genome
fi
echo

# Sorting by coordinates and indexing
echo "Sorting by coordinates and indexing"
echo
# MirBase alignments
echo "Sorting and indexing MirBase alignments"
samtools sort -@ 16 -T $outputdir/"tmp" $outputdir/$sample/$sample"-maturemiRNA-aligned-bowtie1-merged.bam" > $outputdir/$sample/$sample"-maturemiRNA-aligned-bowtie1-CoordSort.bam"
samtools index -@ 16 $outputdir/$sample/$sample"-maturemiRNA-aligned-bowtie1-CoordSort.bam"
rm -rf $outputdir/$sample/$sample"-maturemiRNA-aligned-bowtie1-merged.bam"
#rm -rf $outputdir/$sample/$sample"-maturemiRNA-aligned-bowtie1-fixmate.bam"
echo Done
echo
# Reference genome alignments
echo "Sorting and indexing genome reference alignments for downstream analyses"
samtools sort -@ 16 -T $outputdir/"tmp" $outputdir/$sample/$sample"-genomeaftermiRNA-aligned-bowtie1-merged.bam" > $outputdir/$sample/$sample"-genomeaftermiRNA-aligned-bowtie1-CoordSort.bam"
samtools index -@ 16 $outputdir/$sample/$sample"-genomeaftermiRNA-aligned-bowtie1-CoordSort.bam"
rm -rf $outputdir/$sample/$sample"-genomeaftermiRNA-aligned-bowtie1-merged.bam" 
#rm -rf $outputdir/$i/$i"-genomeaftermiRNA-aligned-bowtie1-fixmate.bam"
echo "Sorting and indexing done"
echo

#echo "Counting step"
#This counting step is for libraries where amplification is made before library preparation as in single cell
# For bulk preparation use umi_tools dedupl
#echo "Couting step without UMI's will not use Umi tools. Continue counting with samtools idxstats"
#umi_tools count --method=unique --per-contig -I $outputdir/$i/$i"-maturemiRNA-aligned-bowtie1-CoordSort.bam" \
#-L $outputdir/$i/$i"_counts-uniquemethod-maturemiRNA.log" -S $outputdir/$i/$i"_counts-finaloutput-uniquemethod-maturemiRNA.txt"
#echo

echo "Deduplicating the aligned maturemiRNA BAM"
# echo "Removal of duplicated read isn't recomended for RNAseq data without UMIs"
umi_tools dedup --method=unique -I $outputdir/$sample/$sample"-maturemiRNA-aligned-bowtie1-CoordSort.bam" \
-S $outputdir/$sample/$sample"_deduplicated-matureMirna-uniquemethod-beststratam1.bam" -L $outputdir/$sample/$sample"-deduplicate-matureMirna-uniquemethod-beststratam1.log"
echo
echo "Indexing the BAM output for finding counts"
#echo "Not needed for data without UMIS"
samtools index -@ 16 $outputdir/$sample/$sample"_deduplicated-matureMirna-uniquemethod-beststratam1.bam"
echo

echo "Deduplicating the aligned genome aligned sortedBAM"
#echo "Removal of duplicated read isn't recomended for RNAseq data without UMIs"
umi_tools dedup --method=unique -I $outputdir/$sample/$sample"-genomeaftermiRNA-aligned-bowtie1-CoordSort.bam" \
-S $outputdir/$sample/$sample"_deduplicated-genomeaftermiRNA-uniquemethod.bam" -L $outputdir/$sample/$sample"-deduplicate-genomeaftermiRNA-uniquemethod.log"
echo
echo "Indexing the genome aligned BAM output for finding counts"
echo "Not needed for data without UMIS"
samtools index -@ 16 $outputdir/$sample/$sample"_deduplicated-genomeaftermiRNA-uniquemethod.bam"
echo

echo "MirBase miRNA counting for sample "${sample}
samtools idxstats -@ 16 $outputdir/$sample/$sample"_deduplicated-matureMirna-uniquemethod-beststratam1.bam" \
| cut -f1,3 - | sed "1s/^/miRNA\t${i}-miRNAcount\n/" - > $outputdir/maturemiRNAcounts/${sample}"-maturemiRNA-counts.txt"
echo Done
echo

echo "Genome aligned miRNA counting for sample "${sample}
echo "Actual filtering sam step only for overlapping microRNA locations"
tagBam -i $outputdir/$sample/$sample"_deduplicated-genomeaftermiRNA-uniquemethod.bam" -files $humanMirTagBamBed \
-names -tag XQ > $outputdir/$sample/$sample"-tagged.bam"
echo Done
echo

echo "Completed analysis for whole run!!"
