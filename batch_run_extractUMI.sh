#!/bin/bash
#SBATCH --job-name=remove_adapters
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=20G


##################
# Create block to find and split bam files on gene name
##################

# create destination dir
#mkdir Allbams-$GENEID

# Begin by creating a list of datasets that we want to collect bams for
find . -maxdepth 1 -iname "*R1_001.adapterTrim.round2.fastq" | awk -F'[/]' '{print $2 }' > umi_datasets.txt

declare -a datasets=( $(cut -b 1- ./umi_datasets.txt) )

my_func() {
        sample_id_R1=$(echo $1 | awk -F'[.]' '{print $1 "." $2 "." $3}')
        sample_id_R2=$(echo $sample_id_R1 | awk -F'[_]' '{print $1 "_" $2 "_" $3 "_" $4 "_" $5 "_" "R2" "_" $7 }')

        python eCLIP_ExtractUMI.py \
		-f $sample_id_R1.fastq \
		-r $sample_id_R2.fastq \
		-l "$(wc -l $sample_id_R1.fastq | awk '{print $1}')" \
		-o ./ \
		-s 100000
}


export -f my_func
parallel my_func ::: ${datasets[@]}
