#!/usr/bin/R
parameters=commandArgs(trailingOnly=TRUE) # get command line arguments

usage_statement="
Usage:
	Rscript find_complementary_flips.R <Flipscan Output> <LD Threshold> <Number SNPs> <Output File>
Options:
	Flipscan Output  file created by the plink --flip-scan command (should have extension '.flipscan', edited with 'sed -i 's/NA *$/NA NA/g' to allow R to read it in)
	LD Threshold	 R^2 threshold to use for LD calculations (NEG_R > this number, 0-1)
	Number SNPs	 threshold for the number of SNPs negatively associated with the SNP in another dataset (NEG > this number)
	Output File	 absolute path to desired output file (will overwrite)
Example: Rscript find_complementary_flips.R 2015-9017_2016-9174_filtered.flipscan 0.5 1 2015-9017_2016-9174_remaining_strandedness_problems.snplist"

if(length(parameters)!=4) {
	stop(paste("All calls to this script require you to specify all parameters.", usage_statement))
}

flipscan_output <- read.table(header = TRUE, stringsAsFactors = FALSE, parameters[1])
ld_threshold <- parameters[2]
number_snps <- parameters[3]
output_file <- parameters[4]

check_if_strand_flip <- function(A1,A2) {
  if((A1=="A" & A2=="T")|(A1=="T"& A2=="A"))
    return(TRUE)
  else if((A1=="C" & A2=="G")|(A1=="G" & A2=="C"))
    return(TRUE)
  else
    return(FALSE)
}

strand_errors <- subset(flipscan_output, (NEG > number_snps) & !is.na(R_NEG) & (R_NEG > ld_threshold) & mapply(check_if_strand_flip, A1, A2))
#write.table(strand_errors, file="full_strand_errors.txt", row.names=FALSE, quote=FALSE)
write.table(strand_errors$SNP, file=output_file, row.names=FALSE, quote=FALSE, col.names=FALSE)  
