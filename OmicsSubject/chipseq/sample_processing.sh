#!/bin/bash
#SBATCH --job-name=grupo4
#SBATCH --export=ALL

## Chip-Seq -- Data processing --> Only for experiments using input as control.

sampledir=$1 # Text file with directory of every sample (experimental design is required prior script execution)
# Important --> Last line in the text file must be the reference genome directory
number=$(wc -l < $sampledir) # Number of samples (including control samples)
sampleSRA=$2 # Text file with SRA accession number (One accesion number for each line)
genomedir=$(tail -n 1 < $sampledir)

## Generate reference genome index
cd $genomedir
bowtie2-build *.fa index # Reference genome must be in fasta format

## Processing start
txtsdir=$(cat $3) # Folder with sampledir and sampleSRA routes
cd "$txtsdir"

## Data analysis
for i in $(seq 1 $(($number-1))); do # All the directories are read except reference genome directory
        samplecd=$(sed -n "${i}p" $sampledir) # Assign working directory of the sample
        samplefastq=$(sed -n "${i}p" $sampleSRA) # Set accesion number

        ## Download of the sequence
        cd $samplecd # Set working directory
        fastq-dump --gzip --split-files $samplefastq # Download sequence
        echo "$samplefastq downloaded correctly --------------------------- sample $i" ## Check control
        echo "-----------------------------------------------------------------------"

        ## Quality control
        fastqc "$samplefastq"_1.fastq.gz
        echo "Quality control of $samplefastq done --------------------------- sample $i" ## Check control
        echo "-----------------------------------------------------------------------"

        ## Read mapping to reference genome
        bowtie2 -x $genomedir/index -U "$samplefastq"_1.fastq.gz -S "$samplefastq"_1.sam ## IMPORTANT --> ASSUMING SINGLE READS
        echo "$samplefastq mapped --------------------------- sample $i" ## Check control
        echo "-----------------------------------------------------------------------"

        ## Generting sorted bam file
        samtools sort -o "$samplefastq"_1.bam "$samplefastq"_1.sam # Bam file generated
        rm "$samplefastq"_1.sam # (optional, highly recommended)
        # rm *.fastq.gz # fastq.gz removal (optional)
        samtools index "$samplefastq"_1.bam # Bam index generation
        bamCoverage -bs 10 --normalizeUsing CPM --bam "$samplefastq"_1.bam -o "$samplefastq"_1.bw
        # rm *.bam # bam removal (optional)
        echo "Bam, bam index and bw files for $samplefastq generated" ## Check control
        echo "-----------------------------------------------------------------------"

        cd "$txtsdir" ## Return to folder directory
done

echo "Job finished"
echo "Next step --> Peaks Calling"
