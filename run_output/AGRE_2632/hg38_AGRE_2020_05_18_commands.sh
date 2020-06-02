#!/usr/bin/env bash

# I changed quite a few things, so here are the parameters for running merge_and_filter.sh in May 2020

project_directory="/home/sjensen/AGRE_2632"
output_name="hg38_AGRE"
sample_info="${project_directory}/AGRE_sample_info.ped" 
datasets_file="${project_directory}/dataset_locations_2020_02_27.txt"
snp_names="/home/sjensen/dbsnp_references/dbsnp_153_SNPs_reformatted_reference.txt" #"${project_directory}/hg19_dbsnp144_ACE_rs_ids.txt" # for now just use file generated for just ACE samples
log_file="${project_directory}/${output_name}_merge_and_filter.log"
build=hg38

# QC parameters
# Following Jill Haney's PRS preparation steps (https://docs.google.com/presentation/d/1g4ugFonTy_hWaZAW8JpwHTP-QHtnCxZnhyjoRLdhN04/edit?usp=sharing) 
# and the Michigan Imputation Server guidelines (https://imputationserver.sph.umich.edu/index.html#!pages/help).

# decided to remove QC filters that will limit the number of variants used for imputation
maf=0.01
geno=0.1
hwe=1e-40
mind=0.2

# no MAF or HWE filtering
#/home/sjensen/merge_and_filter.sh --output_directory $project_directory --output_name ${output_name}_geno_0.1_mind_0.2 --sample_info $sample_info --datasets_file $datasets_file --snp_names $snp_names --build $build --geno $geno --mind $mind | tee -a $log_file

# lower MAF threshold
/home/sjensen/merge_and_filter.sh --output_directory $project_directory --output_name ${output_name}_maf_${maf}_geno_${geno}_hwe_${hwe}_mind_${mind} --sample_info $sample_info --datasets_file $datasets_file --snp_names $snp_names --build $build --maf $maf --geno $geno --hwe $hwe --mind $mind | tee -a $log_file

# gzip -c ACE2.vcf > ACE2.vcf.gz
# move gzipped file to Hoffman2 /u/home/s/samlj/project-geschwind/ACE/ACE2.vcf.gz
