#!/usr/bin/R
parameters=commandArgs(trailingOnly=TRUE) # get command line arguments

usage_statement="Usage:
        Rscript rename_variants.R <PLINK dataset SNPs> <Matching dbSNP SNPs> <Output prefix>
Options:
        PLINK dataset SNPs	  Tab delimited file describing SNPs from a PLINK dataset with format 'chromosome position rsID A1 A2' for each variant (no header).
	Matching dbSNP SNPs	  Tab delimited file with dbSNP variants that match the PLINK dataset SNPs by position, with the same format as above.
	Output prefix		  Prefix of files matched variants of each type will be written to.
Example: Rscript rename_variants.R 2017-9154_snps.txt 2017-9154_hg38_matching_snps.txt 2017-9154_hg38_names_and_alleles"

if(length(parameters)!=3) {
        stop(paste("All calls to this script require you to specify all parameters.", usage_statement))
}

library("dplyr")
library("tidyr")

complements <- c("A"="T", "T"="A", "G"="C", "C"="G")

current_snps <- read.table(stringsAsFactors = FALSE, parameters[1]) %>%
	rename(chromosome=V1, position=V2, old_rs_id=V3, A1=V4, A2=V5)

matching_dbsnp <- read.table(stringsAsFactors = FALSE, parameters[2]) %>%
	rename(chromosome=V1, position=V2, dbsnp_rs_id=V3, reference=V4, alternates=V5) %>%
	mutate(alternates=strsplit(alternates,split=",")) %>%
	inner_join(current_snps) %>%
	mutate(complement_A1=complements[A1], complement_A2=complements[A2])


output_prefix=parameters[3]

true_matches <- subset(matching_dbsnp, reference==A1 & A2 %in% alternates) %>% subset(!duplicated(old_rs_id)) # If variant has more than one rs_id, use first.
flipped_matches <- subset(matching_dbsnp, reference==A2 & A1 %in% alternates) %>% subset(!duplicated(old_rs_id))
complement_matches <- subset(matching_dbsnp, complement_A1==reference & complement_A2 %in% alternates) %>% subset(!duplicated(old_rs_id))
complement_flipped_matches <- subset(matching_dbsnp, complement_A2==reference & complement_A1 %in% alternates) %>% subset(!duplicated(old_rs_id))

true_matches_by_snp <- data.frame(table(true_matches$old_rs_id)) %>% rename(old_rs_id=Var1, true_matches=Freq)
flipped_matches_by_snp <- data.frame(table(flipped_matches$old_rs_id)) %>% rename(old_rs_id=Var1, flipped_matches=Freq)
complement_matches_by_snp <- data.frame(table(complement_matches$old_rs_id)) %>% rename(old_rs_id=Var1, complement_matches=Freq)
complement_flipped_matches_by_snp <- data.frame(table(complement_flipped_matches$old_rs_id)) %>% rename(old_rs_id=Var1, complement_flipped_matches=Freq)


current_snps_type <- mutate(current_snps, 
                            true_match=(old_rs_id %in% true_matches$old_rs_id), 
                            flipped_match=(old_rs_id %in% flipped_matches$old_rs_id),
                            complement_match=(old_rs_id %in% complement_matches$old_rs_id),
                            complement_flipped_match=(old_rs_id %in% complement_flipped_matches$old_rs_id)) %>%
  full_join(true_matches_by_snp) %>%
  full_join(flipped_matches_by_snp) %>%
  full_join(complement_matches_by_snp) %>%
  full_join(complement_flipped_matches_by_snp) %>%
  mutate(true_matches=replace_na(true_matches, 0),
         flipped_matches=replace_na(flipped_matches, 0),
         complement_matches=replace_na(complement_matches, 0),
         complement_flipped_matches=replace_na(complement_flipped_matches, 0),
         total_matches=true_matches+flipped_matches+complement_matches+complement_flipped_matches,
	 no_match=(total_matches==0))


# OUTPUT FILE WITH ALL INFO FOR LATER GRAPHING

all_info_name <- paste0(output_prefix, "_info.txt")
write.table(current_snps_type, all_info_name, row.names=FALSE)

# OUTPUT FILE WITH COUNTS FOR QUICK CHECKS

current_snps_by_type <- select(current_snps_type,true_match,flipped_match,complement_match,complement_flipped_match,no_match)
all_counts <- data.frame(table(current_snps_by_type)) %>% subset(Freq != 0)
all_counts_name <- paste0(output_prefix, "_counts.txt")
write.table(all_counts, all_counts_name, row.names=FALSE, quote=FALSE)

# OUTPUT FILE WITH ALL VARIANTS WITH NO MATCHES (TO BE REMOVED)

no_matches <- subset(current_snps_type, total_matches==0) %>%
	select(old_rs_id)
no_match_output_name <- paste0(output_prefix, "_no_match.txt")
write.table(no_matches, no_match_output_name, row.names=FALSE, col.names=FALSE, quote=FALSE)

# OUTPUT FILE WITH ALL TRUE MATCHES (TO BE RENAMED)

true_matches_output <- select(true_matches, old_rs_id, dbsnp_rs_id)
true_match_output_name <- paste0(output_prefix, "_true_matches.txt")
write.table(true_matches_output, true_match_output_name, row.names=FALSE, col.names=FALSE, quote=FALSE)

# OUTPUT FILE WITH REMAINING FLIPPED MATCHES

flipped_matches_output <- subset(flipped_matches, !(old_rs_id %in% true_matches$old_rs_id)) %>%
	select(old_rs_id, dbsnp_rs_id, reference)
flipped_matches_output_name <- paste0(output_prefix, "_flipped_matches.txt")
write.table(flipped_matches_output, flipped_matches_output_name, row.names=FALSE, col.names=FALSE, quote=FALSE)

# OUTPUT REMAINING COMPLEMENT MATCHES

complement_matches_output <- subset(complement_matches, !(old_rs_id %in% true_matches$old_rs_id) & !(old_rs_id %in% flipped_matches$old_rs_id)) %>%
	select(old_rs_id, dbsnp_rs_id)
complement_matches_output_name <- paste0(output_prefix, "_complement_matches.txt")
write.table(complement_matches_output, complement_matches_output_name, row.names=FALSE, col.names=FALSE, quote=FALSE)

# OUTPUT REMAINING COMPLEMENT FLIPPED MATCHES

complement_flipped_matches_output <- subset(complement_flipped_matches, !(old_rs_id %in% true_matches$old_rs_id) & !(old_rs_id %in% flipped_matches$old_rs_id) & !(old_rs_id %in% complement_matches$old_rs_id)) %>%
	select(old_rs_id, dbsnp_rs_id, reference)
complement_flipped_matches_output_name <- paste0(output_prefix, "_complement_flipped_matches.txt")
write.table(complement_flipped_matches_output, complement_flipped_matches_output_name, row.names=FALSE, col.names=FALSE, quote=FALSE)

