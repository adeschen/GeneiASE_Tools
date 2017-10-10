#!/usr/bin/env Rscript

library("optparse")

##############################################################
## Arguments parsing
##############################################################
option_list = list(
    make_option(c("-f", "--file"), type="character", action="store", default=NULL, 
                help="Filtered SNPs file without header (output of GATK VariantFilter) (header lines starting with # will be discarded).\n\t\t The 1st and 2nd columns should be for the chromosome and the position.\n\t\t The 3th and 4th columns should be for the reference and alternative allele.\n\t\t The 7th column should contain the filtering results.", metavar="character"),
    make_option(c("-v", "--vep"), type="character", action="store", default=NULL, 
                help="Vep annotated file,\n\t\t The Ensembl Gene ID must be in the 4th column.", metavar="character"),
    make_option(c("-m", "--matrix"), type="character", action="store", default=NULL, 
                help="Depth of alignment file", metavar="character"),
    make_option(c("-n", "--nomatrix"), type="character", action="store", default=NULL, 
                help="Rejected depth of alignment file for specific filtering conditions.", metavar="character"),
    make_option(c("-d", "--depth"), type="integer", action="store", default=10, 
                help="Minimum depth of alignment (reference + alternative allele)", metavar="number"),
    make_option(c("-o", "--out"), type="character", action="store", default="out.txt", 
                help="output file name [default= %default]", metavar="character")
); 

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$file) || is.null(opt$vep) || is.null(opt$matrix) || is.null(opt$nomatrix)){
    print_help(opt_parser)
    stop("At least 5 arguments must be supplied.\n", call.=FALSE)
}


##############################################################
## Load VEP annotation file
##############################################################

## File with gene name (only keep when Ensembl Gene ID present)
## Chr:Position must be in the 2nd column
## Ensembl ID must be in the 4th column
vep <- read.table(opt$vep, header=F, stringsAsFactors = FALSE)
vep <- subset(vep, V4 != "-")


##############################################################
## Load Depth files and extract filtered depth of alignement
##############################################################

## File with depth information (all reads)
all <- read.table(opt$matrix, header=TRUE, stringsAsFactors = FALSE)

## File with depth rejected by filtering
noGood <- read.table(opt$nomatrix, header=TRUE, stringsAsFactors = FALSE)

## Generate table with depth related to filtering conditions
## Real depth for filtering = all - depth rejected by filtering
## all - noGood
goodDepth <- noGood
goodDepth[, 3:6] <- all[,3:6] - noGood[,3:6]
goodDepth$ID <- paste0(goodDepth[,1], ":", goodDepth[,2])
rm(all)
rm(noGood)


##############################################################
## Load SNPs file
##############################################################

## Filter SNPs file
allSNPs <-  read.table(opt$file,  header=F, stringsAsFactors = FALSE)
## Only PASS SNPs are retained
passSNPs <- subset(allSNPs, V7 == "PASS")
rm(allSNPs)
## Remove indels (everything that is more than 1 pair base)
passSNPs <- subset(passSNPs, nchar(V4) == 1)
passSNPs <- subset(passSNPs, nchar(V5) == 1)
## Remove extra columns
passSNPs <- passSNPs[, 1:5]
## Create unique NAME : chr + ":" + position
passSNPs[,6] <- paste0(passSNPs[,1], ":", passSNPs[,2])
colnames(passSNPs) <- c("CHRM", "POS", "ID", "REF", "ALT", "NAME")


##############################################################
## Add depth informations to the SNP information using CHR+POSITION AS KEY
##############################################################

## Combine the information from SNPs and depth
newSNPs <- merge(passSNPs, goodDepth, by.x = "NAME", by.y = "ID")
newSNPs$ALT_COUNT <- rep(-1, length(newSNPs$NAME))
newSNPs$REF_COUNT <- rep(-1, length(newSNPs$NAME))

## Selected depth for the reference and alternative allele
for (letter in c("A", "C", "G", "T")) {
    alt_positions <- newSNPs$ALT == letter
    newSNPs$ALT_COUNT[alt_positions] <- newSNPs[, letter][alt_positions]
    ref_positions <- newSNPs$REF == letter
    newSNPs$REF_COUNT[ref_positions] <- newSNPs[, letter][ref_positions]
}

## Retained only needed columns
newSNPs <- newSNPs[,c("NAME", "ALT_COUNT", "REF_COUNT")]


##############################################################
## Add depth informations to the SNP information using CHR+POSITION AS KEY
##############################################################

## Combine the depth information with the gene information
newFinal <- merge(newSNPs, vep, by.x = "NAME", by.y="V2")
newFinal <- newFinal[, c("NAME", "ALT_COUNT", "REF_COUNT", "V4")]
newFinal <- newFinal[!duplicated(newFinal), ]
newFinal <- subset(newFinal, ALT_COUNT + REF_COUNT >= opt$depth)

final <- newFinal[, c("V4", "NAME", "ALT_COUNT", "REF_COUNT")]
colnames(final) <- c("gene", "snp.id", "alt.dp", "ref.dp")

## Save result in output file
write.table(final, file = opt$out, quote=FALSE, sep="\t", 
                row.names = FALSE, col.names = TRUE)

