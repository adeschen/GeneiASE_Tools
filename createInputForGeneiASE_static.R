#!/usr/bin/env Rscript

library("optparse")

## Arguments parsing
option_list = list(
    make_option(c("-f", "--file"), type="character", action="store_true", default=NULL, 
                help="Filtered SNPs file (output of GATK VariantFilter)", metavar="character"),
    make_option(c("-v", "--vep"), type="character", action="store_true", default=NULL, 
                help="Vep annotated file", metavar="character"),
    make_option(c("-m", "--matrix"), type="character", action="store_true", default=NULL, 
                help="Depth of alignment file", metavar="character"),
    make_option(c("-n", "--nomatrix"), type="character", action="store_true", default=NULL, 
                help="Rejected depth of alignment file for specific filter", metavar="character"),
    make_option(c("-d", "--depth"), type="integer", action="store_true", default=10, 
                help="Minimum depth of alignment (reference + alternative allele)", metavar="number"),
    make_option(c("-o", "--out"), type="character", action="store_true", default="out.txt", 
                help="output file name [default= %default]", metavar="character")
); 

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$file) || is.null(opt$vep) || is.null(opt$matrix) || is.null(opt$nomatrix)){
    print_help(opt_parser)
    stop("At least 5 arguments must be supplied.\n", call.=FALSE)
}

## File with gene name (only keep when gene name present)
vep <- read.table(opt$vep, header=F, stringsAsFactors = FALSE)
vep <- subset(vep, V4 != "-")

## File with depth information (all reads)
all <- read.table(opt$matrix, header=TRUE, stringsAsFactors = FALSE)

## File with depth rejected by filtering
## Real depth for filtering = all - depth rejectect by filtering
## all - noGood
noGood <- read.table(opt$nomatrix, header=TRUE, stringsAsFactors = FALSE)

## Generate table with depth related to filtering conditions
goodDepth <- noGood
goodDepth[, 3:6] <- all[,3:6] - noGood[,3:6]
goodDepth[,7] <- paste0(goodDepth[,1], ":", goodDepth[,2])

## Filter SNPs file. 
allSNPs <-  read.table(opt$file,  header=F, stringsAsFactors = FALSE)
## Only PASS SNPs are retained
passSNPs <- subset(allSNPs, V7 == "PASS")
## Remove indels
passSNPs <- subset(passSNPs, nchar(V4) == 1)
passSNPs <- subset(passSNPs, nchar(V5) == 1)
## Remove extra columns
passSNPs <- passSNPs[, 1:5]
## Create unique NAME : chr + position
passSNPs[,6] <- paste0(passSNPs[,1], ":", passSNPs[,2])
colnames(passSNPs) <- c("CHRM", "POS", "ID", "REF", "ALT", "NAME")


## Combine the information from SNPs and depth
newSNPs <- merge(passSNPs, goodDepth, by.x = "NAME", by.y = "V7")
newSNPs$ALT_COUNT <- rep(-1, length(newSNPs$NAME))
newSNPs$REF_COUNT <- rep(-1, length(newSNPs$NAME))

## Selected depth for the reference and alternative allele
for (i in seq_len(nrow(newSNPs))) {
    newSNPs$ALT_COUNT[i] <- newSNPs[i, newSNPs$ALT[i]]
    newSNPs$REF_COUNT[i] <- newSNPs[i, newSNPs$REF[i]]
}

## Rename columns
final <- newSNPs[,c("NAME", "ALT_COUNT", "REF_COUNT")]

## Combine the depth information with the gene information
newFinal <- merge(final, vep, by.x = "NAME", by.y="V2")
newFinal <- newFinal[, c("NAME", "ALT_COUNT", "REF_COUNT", "V4")]
newFinal <- newFinal[!duplicated(newFinal), ]
newFinal <- subset(newFinal, ALT_COUNT + REF_COUNT >= opt$depth)

final <- newFinal[, c("V4", "NAME", "ALT_COUNT", "REF_COUNT")]
colnames(final) <- c("gene", "snp.id", "alt.dp", "ref.dp")

## Save result in output file
write.table(final, file = opt$out, quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)

