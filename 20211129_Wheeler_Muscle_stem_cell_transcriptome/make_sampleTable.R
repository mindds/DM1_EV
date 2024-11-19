#!/usr/bin/R

# this program is used to run merge lanes by R1 & R2
# 03/10/2022 Zhaozhi Li

# usage: Rscript make_sampleTable.R -f fastq_dir -o sampleTable_name
#   fastaq_dir: path to the folder that has all fastq files
#   sampleTable_name: path to sampleTable.csv

# reference:
# https://www.biostars.org/p/184085/

# logs
# 03/10/2022 
#   Rscript make_sampleTable.R -f /space/mindds/1/projects/ThurmanRNASeq/data/merged/20211129_Wheeler_Muscle_stem_cell_transcriptome
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# set up
suppressMessages({
  require(optparse)
  require(tidyr)
  require(data.table)
})

# get options
option_list = list(
  make_option(c("-f", "--fastqdir"), action="store", type="character", default=NULL, 
              help="path to merged fastq files", metavar="character"),
  make_option(c("-o", "--out"), action="store", type="character", default="sampleTable.csv", 
              help="filename of sampleTable", metavar="character")
)

opt = parse_args(OptionParser(option_list=option_list))
if (is.null(opt$fastqdir)){
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}


#-------------------------------------------------------------------------------
# generate sampleTable.csv

tlb = data.table(fastqName = list.files(opt$fastqdir))
# sampleID from fastq.gz file names
tlb[, sampleID := gsub("_R[12]_001.fastq.gz", "", fastqName)]
# R1/R2 from fastq.gz file names
tlb[, read := ifelse(grepl("_R1_001.fastq", fastqName), "R1", 
              ifelse(grepl("_R2_001.fastq", fastqName), "R2", ""))]
if(any(is.na(tlb$read))){
  stop("can't distinguish R1/R2 from fastq.gz\n")
}
tlb = tlb[, .(sampleID, read, fastqName)]

# output
fwrite(tlb, file = opt$out)
