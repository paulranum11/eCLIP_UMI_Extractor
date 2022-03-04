#!/bin/bash
#SBATCH --job-name=remove_adapters
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G


##################
# Create block to find and split bam files on gene name
##################

# create destination dir
#mkdir Allbams-$GENEID

# Begin by creating a list of datasets that we want to collect bams for
find . -maxdepth 1 -iname "*eCLIP*R1*" | awk -F'[/]' '{print $2 }' > datasets.txt

declare -a datasets=( $(cut -b 1- ./datasets.txt) )

my_func() {
        sample_id_R1=$(echo $1 | awk -F'[.]' '{print $1 }')
        sample_id_R2=$(echo $sample_id_R1 | awk -F'[_]' '{print $1 "_" $2 "_" $3 "_" $4 "_" $5 "_" "R2" "_" $7 }')

        python eCLIP_ExtractUMI.py \
		-f $sample_id_R1.fastq \
		-r $sample_id_R2.fastq \
		-l "$(wc -l $sample_id_R1.fastq | awk '{print $1}')" \
		-o ./ \
		-s 5
}


export -f my_func
parallel my_func ::: ${datasets[@]}
