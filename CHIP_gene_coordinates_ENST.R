# this script is to get the ENST ID and coordinates for genes which have CHIP
setwd("D:/2022-04_Mutect/mutect_CHIP")

library(stringi)   

CHIP_mutations <- read.xlsx("D:/2022-04_Mutect/mutect_CHIP/Supplementary Table 2. Leukemogenic driver mutations queried in this study.xlsx")

# find ENST ID
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("biomaRt")
library(biomaRt)
#detach(package:biomaRt)
ensembl <-  useMart("ensembl", dataset="hsapiens_gene_ensembl")
grch37 <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", dataset="hsapiens_gene_ensembl")

## getting genomic coordiante for CHIP genes
attributes <- c("ensembl_transcript_id","start_position","end_position","strand","hgnc_symbol","chromosome_name", "refseq_mrna")
CHIP_gene_coordinates <- getBM(attributes=attributes, 
                               filters = "refseq_mrna", values = CHIP_mutations$Accession, mart= grch37)
write.table(CHIP_gene_coordinates, file = "CHIP_gene_coordinates_ENST.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)

detach(package:biomaRt)
