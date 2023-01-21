#!/bin/bash
#SBATCH --job-name=index
#SBATCH --export=ALL
#SBATCH --output=genome_index

## Chip-Seq -- Data processing --> Only for experiments using input as control.

sampledir=$1 # Text file with directory of every sample (experimental design is required prior script execution)
# Important --> Last line in the text file must be the reference genome directory
number=$(wc -l < $sampledir) # Number of samples (including control samples)
sampleSRA=$2 # Text file with SRA accession number (One accesion number for each line)

## Processing start
for i in $(seq 1 $(($number-1))); do # All the directories are read except reference genome directory
        samplecd=$(sed -n "${i}p" $sampledir) # Assign working directory of the sample
        samplefastq=$(sed -n "${i}p" $sampleSRA) # Set accesion number

## Download of the sequence
        cd $samplecd # Set working directory
        fastq-dump --gzip --split-files $samplefastq # Download sequence
        echo "$samplefastq downloaded correctly" ## Check control

## Quality control
        fastqc $samplefastq.fastq.gz
        echo "Quality control of $sampleSRA done" ## Check control

## Read mapping to reference genome
        genomedir=$(tail -n 1 < $sampledir)

        cd $genomedir
        bowtie2-build ?.fa index # Reference genome must be in fasta format

        cd $samplecd
        bowtie2 -x $genomedir/index -U $samplefastq.fastq.gz -S $samplefastq.sam ## IMPORTANT --> ASSUMING SINGLE READS

        echo "$sampleSRA mapped" ## Check control

## Generting sorted bam file
        samtools sort -o $samplefastq.bam $samplefastq.sam # Bam file generated
        rm $samplefastq.sam # (optional, highly recommended)
        # rm *.fastq.gz # fastq.gz removal (optional)
        samtools index $samplefastq.bam # Bam index generation
        bamCoverage -bs 10 --normalizeUsing CPM --bam $samplefastq.bam -o $samplefastq.bw
        # rm *.bam # bam removal (optional)
        echo "Bam, bam index and bw files for $samplefastq generated" ## Check control
done

echo "Job finished"
echo "Next part --> Peaks Calling"
