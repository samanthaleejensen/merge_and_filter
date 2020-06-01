#!/usr/bin/env bash

awk '{print $1}' dbsnp144_plink_mapped_variants_with_alleles.txt | uniq -c | while read num dupe; do [[ $num > 1 ]] && grep -- "^${dupe}\t" dbsnp144_plink_mapped_variants_with_alleles.txt; done
