#!/bin/bash
#$ -cwd
#$ -N get_dbsnp_153_variant_info
#$ -j y
#$ -l h_data=8G,h_rt=24:00:00
#$ -m ea

### GET ONLY QC-PASSED VARIANT POSITIONS

# run QC steps on 2017-9154_2016-9174-3_2016-9174_2015-9017_2013-111A_2010-002_CHOP_2012-092_final using just_qc.sh script
# move results from plink_files to this directory
# delete all but bim
# rename bim to ACE2_variants.bim
# get just locations of variants from second column of bim file
#awk '{print $2}' ACE2_variants.bim > ACE2_snps.txt


### FORMAT DBSNP VERSION 153 (DOWNLOADED FROM NCBI at https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.38.gz)

# 1. Get only Primary Assembly SNVs
dbsnp_vcf=/u/home/s/samlj/project-geschwind/ACE/genotypes/orion_code/get_dbsnp_names/dbsnp_153/GCF_000001405.38.gz
dbsnp_snvs=dbsnp_153_SNVs_reformatted.txt #dbsnp_153_SNVs.vcf

# SNVs will have VC=SNV in the INFO field
# Variants that are part of the primary assembly (AKA chr 1-22, X, Y, MT) all have "NC_" at the beginning of their refseq sequences, which is what's in the chromosome position of these VCFs for some reason.

#zcat $dbsnp_vcf | grep "VC=SNV" | grep "NC_" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > $dbsnp_snvs

# 2. Reformat weird chromosome IDs 

# Since NCBI is dumb, the chromosome column is not the chromosome number but the RefSeq access number, so I have to convert them to chromosomes. I downloaded the reference from https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_assembly_report.txt and removed all lines except those for the primary assembly and mitochondrial DNA. I then kept only the column with the chromosome name and RefSeq accession number, saving that output to GRCh38_primary_assembly_reference.txt.

primary_assembly_reference="GRCh38_primary_assembly_reference.txt"

#grep "NC_" GCF_000001405.38_GRCh38.p12_assembly_report.txt | awk '{print $7"\t"$1}' > $primary_assembly_reference

all_dbsnp_153_SNPs="dbsnp_153_SNPs_reformatted_reference.txt"
while read -r -a chromosome_reference; do
	refseq=${chromosome_reference[0]}
	chromosome=${chromosome_reference[1]}
	echo "Converting $refseq to $chromosome in dbsnp153 documentation."
	grep $refseq $dbsnp_snvs | sed "s/$refseq/$chromosome/g" >> $all_dbsnp_153_SNPs
done < $primary_assembly_reference

#sed '/_/d' $dbsnp_snvs
#awk '{print $1":"$3"\t"
#sed "1D" dbsnp144_all_snps.txt | awk '{print $1":"$3"\t"$4}'| sed 's/chr//g' > dbsnp144_plink_snp_resource.txt
#sed '/_/d' dbsnp144_all_snps.txt > dbsnp144_mapped_variants.txt
#awk ' $6 != "" ' dbsnp144_mapped_variants.txt > dbsnp144_mapped_variants_with_alleles.txt
#sed "1D" dbsnp144_mapped_variants_with_alleles.txt | awk '{print $1":"$3"\t"$4"\t"$5"\t"$6"\t"$7}'| sed 's/chr//g' > dbsnp144_plink_mapped_variants_with_alleles.txt


### GET NAMES THAT MATCH EACH QC-PASSED VARIANT POSITION

#grep -F -f ACE2_snps.txt dbsnp144_plink_mapped_variants_with_alleles.txt > ACE_variant_names.grep.txt

### USE hg19_dbSNP_parsing.R SCRIPT ON LAPTOP TO DETERMINE WHAT EACH POSITION SHOULD BE NAMED
