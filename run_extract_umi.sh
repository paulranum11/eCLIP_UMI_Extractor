#!/bin/bash
#SBATCH --job-name=remove_adapters
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G

#python eCLIP_ExtractUMI.py -f ./BA9_PTBP2_eCLIP_IP1_S10_R1_001_head.fastq -r ./BA9_PTBP2_eCLIP_IP1_S10_R2_001_head.fastq -l "$(wc -l ./BA9_PTBP2_eCLIP_IP1_S10_R1_001_head.fastq | awk '{print $1}')" -o ./ -b 1000

#linecount="$(wc -l ./head_20_R1.fastq | awk '{print $1}')"
#acceptable_bin="$(((linecount - 1) / 4))"

python eCLIP_ExtractUMI.py -f BA9_PTBP2_eCLIP_IP1_S10_R1_001_head.fastq -r BA9_PTBP2_eCLIP_IP1_S10_R2_001_head.fastq -l "$(wc -l ./BA9_PTBP2_eCLIP_IP1_S10_R1_001_head.fastq | awk '{print $1}')" -o ./ -s 5
