#!/usr/bin/env Rscript

library("optparse")

## Arguments parsing
option_list = list(
    make_option(c("-i", "--input"), type="character", action="store", default = NULL, 
                help="Input file with header and Ensembl Gene Name in first column", metavar="character"),
    make_option(c("-b", "--biomart"), type="character", action="store", default = NULL, 
                help="Biomart file with header, Ensembl Gene Name in first column and Gene Name in third column", metavar="character"),
    make_option(c("-o", "--out"), type="character", action="store", default="out.txt", 
                help="output file [default= %default]", metavar="character")
); 

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$biomart)) {
    print_help(opt_parser)
    stop("At least 2 arguments must be supplied.\n", call.=FALSE)
}

## File with Ensembl Gene Name
input <- read.table(opt$input, header = T, sep="\t", stringsAsFactors = FALSE)

## Biomart file with Ensemble Gene Name and Gene Name
biomart <- read.table(opt$biomart, header = T, sep="\t", stringsAsFactors = FALSE)
biomart <- biomart[, c(1,3)]

## Merge result with biomart data
result <- merge(biomart, input, by.x = 1, by.y = 1, all.x = FALSE, all.y = FALSE)
result <- result[!duplicated(result), ]

## Save merged data
write.table(result, opt$out, quote = FALSE, sep="\t", row.names = FALSE, col.names = TRUE)