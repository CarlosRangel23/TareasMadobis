#!/bin/bash
#SBATCH --job-name=grupo4
#SBATCH --export=ALL

## Chip-Seq -- Data processing --> Only for experiments using input as control.

## Variables used during the script
# Important --> In the experimental design, all the inputs folders must be named as input.
sampledir=$1 # Text file with directory of every sample (experimental design is required prior script execution)
# Important --> In sampledir file, inputs files must be the last ones
number=$(wc -l < $sampledir) # All samples, including inputs
input=$(grep "input" $sampledir | wc -l) # Only inputs route
Ninput=$(grep "input" $sampledir | wc -l) # Number of inputs done
sample=$(grep -v "input" $sampledir) # Only samples route
Nsamples=$(grep -v "input" $sampledir | wc -l) # Number of samples

otherdir=$2 # Text file where: first line: route to reference genome; second line: ssodir; third line: route to results folder
genomedir=$(head -n 1 < $otherdir) # Reference genome directory --> Reference genome must be in fasta format
resultsdir=$(tail -n 1 < $otherdir) # Results directory
ssodir=$(head -n 2 < $otherdir | tail -n 1) # Route to folder with sampledir, otherdir and sampleSRA routes

sampleSRA=$3 # Text file with SRA accession number (One accesion number for each line)

echo " "
echo "#############################"
echo "### STARTING THE ANALYSIS ###"
echo "#############################"
echo " "

## Generate reference genome index
cd $genomedir
bowtie2-build *.fa index # Reference genome must be in fasta format

## Data analysis
cd $ssodir
for i in $(seq 1 $number); do # All the directories are read
        samplecd=$(sed -n "${i}p" $sampledir) # Assign working directory of the sample
        samplefastq=$(sed -n "${i}p" $sampleSRA) # Set accesion number

        ## Download of the sequence
        cd $samplecd # Set working directory
        fastq-dump --gzip --split-files $samplefastq # Download sequence
        echo "-----------------------------------------------------------------------"
        echo "$samplefastq downloaded correctly --------------------------- sample $i" ## Check control
        echo "-----------------------------------------------------------------------"

        ## Quality control
        fastqc "$samplefastq"_1.fastq.gz
        echo "-----------------------------------------------------------------------"
        echo "Quality control of $samplefastq done --------------------------- sample $i" ## Check control
        echo "-----------------------------------------------------------------------"

        ## Read mapping to reference genome
        bowtie2 -x $genomedir/index -U "$samplefastq"_1.fastq.gz -S "$samplefastq"_1.sam ## IMPORTANT: ASSUMING SINGLE READS
        echo "-----------------------------------------------------------------------"
        echo "$samplefastq mapped --------------------------- sample $i" ## Check control
        echo "-----------------------------------------------------------------------"

        ## Generating sorted bam file
        samtools sort -o "$samplefastq"_1.bam "$samplefastq"_1.sam # Bam file generated
        rm "$samplefastq"_1.sam # (optional, highly recommended)
        rm *.fastq.gz # fastq.gz removal (optional)
        samtools index "$samplefastq"_1.bam # Bam index generation
        bamCoverage -bs 10 --normalizeUsing CPM --bam "$samplefastq"_1.bam -o "$samplefastq"_1.bw
        # rm *.bam # bam removal (optional)
        echo "-----------------------------------------------------------------------"
        echo "Bam, bam index and bw files for $samplefastq generated" ## Check control
        echo "-----------------------------------------------------------------------"

        cd $ssodir ## Return to folder directory
done

echo " "
echo "###############################"
echo "#### STARTING PEAK CALLING ####"
echo "###############################"
echo " "

## Data processing only with one input

if [ $Ninput -eq 1 ]; then
        onlyinput=$(head -n 1 < $input)
        for i in $(seq 1 $Nsamples); do # All the chip samples are read
                samplecd=$(grep -v "input" $sampledir | sed -n ${i}p) # Assign directory of the bam chip sample file
                samplefastq=$(sed -n "${i}p" $sampleSRA) # Set accesion number

                cd $resultsdir
                macs2 callpeak -t $samplecd/*.bam -c $onlyinput/*.bam -f BAM --outdir . -n $samplefastq

                echo " "
                echo "-----------------------------------------------------------------------"
                echo "PEAK CALLING FOR $samplefastq DONE" ## Check control
                echo "-----------------------------------------------------------------------"
                echo " "
                cd $ssodir
        done

## Data processing with more than one input
## Here we assume that chip01??s control is input01, for chip02 input02, and so on.
## That is to say, Number of chip samples must be equal to Number of control inputs

else
        if [ $Ninput -ne $Nsamples ]; then
                echo "Fatal error, number of control samples and number of chip samples differ"
                exit
        fi

        for i in $(seq 1 $Nsamples); do # All the chips samples are read
                samplecd=$(grep -v "input" $sampledir | sed -n ${i}p) # Assign directory of the bam chip sample file
                samplefastq=$(sed -n "${i}p" $sampleSRA) # Set accesion number
                inputcd=$(grep "input" $sampledir | sed -n ${i}p) # Assign directory of the bam control input file

                cd $resultsdir
                macs2 callpeak -t $samplecd/*.bam -c $inputcd/*.bam -f BAM --outdir . -n $samplefastq

                echo " "
                echo "-----------------------------------------------------------------------"
                echo "PEAK CALLING FOR $samplefastq DONE" ## Check control
                echo "-----------------------------------------------------------------------"
                echo " "
                cd $ssodir
        done

fi

echo " "
echo "######################"
echo "#### JOB FINISHED ####"
echo "######################"
echo " "
