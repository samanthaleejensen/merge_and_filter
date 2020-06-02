#!/bin/bash
qsub -cwd -N join_ace -e join_ace.log -o join_ace.log -m ae -M samleejensen@gmail.com join_two_files.sh ACE2_snps.txt dbsnp144_plink_mapped_variants_with_alleles.txt ACE_variant_names.txt

qsub -cwd -N join_ace_opposite -e join_ace_opposite.log -o join_ace_opposite.log -m ae -M samleejensen@gmail.com join_two_files.sh dbsnp144_plink_mapped_variants_with_alleles.txt ACE2_snps.txt ACE_variant_names_opposite.txt
 
