setwd("D:/2022-04_Mutect/mutect_CHIP")

library(tidyverse)
library(tidyr)
library(vcfR)
library(broom)
library(openxlsx)
library(data.table)
library(dbplyr)
library(stringi)

readVcf <- function(vcf_file) {
  vcf <- read.vcfR(vcf_file)
  #vcf <- read.vcfR("IPF0017.194167/IPF0017.194167-filtered.ann.vcf.gz")
  v <- vcfR2tidy(vcf, format_fields = c("GT", "AD", "AF", "DP"))
  annd <- v$meta %>% filter(ID == "ANN") %>% select(Description)
  annd <- unlist(str_split(gsub("'", "", gsub(" ","" ,gsub("Functional annotations: ", "", annd$Description))), "\\|"))
  
  ANN <- 	v$fix %>%  
    unite("VariantID", c(CHROM, POS, REF, ALT), sep = "-", remove = FALSE) %>% 
    unite("key", c(ChromKey, POS), sep = "-", remove = FALSE) %>% 
    select("key", "VariantID", "FILTER", "ANN") %>%
    separate_rows(ANN, sep = ",") %>%
    separate(ANN, annd, sep = "\\|", convert = TRUE) %>%
    filter(Feature_ID %in% CHIP_mutations$ensembl_transcript_id) %>% 
    mutate(new.HGVS.p = gsub("p.", "", stri_replace_all_fixed(HGVS.p,
                                                              pattern = paste0(amino_acid_abbr$Abbreviation),
                                                              replacement = paste0(amino_acid_abbr$`Single letter abbreviation`),
                                                              vectorize_all = FALSE)))
  
  gt <- v$gt %>% unite("key", ChromKey, POS, sep = "-", remove = FALSE)
  
  variants <- ANN %>% left_join(gt, by = "key") %>% 
    select(Indiv, everything()) %>% select(-contains(c(".y","Key")))
  
  
  return(variants)
}

amino_acid_abbr <- fread("D:/2022-04_Mutect/amino_acid_abbreviation.txt")
variants.lists <- lapply(list.files(pattern = "-filtered.ann.no-shift.vcf.gz$",recursive = TRUE), function(x) readVcf(x))
variants.list <- do.call(rbind, variants.lists)
variants.list[variants.list == ""] <- " "


GENE_NAME <- "ASXL1"
variants.list %>% 
  filter(str_detect(Annotation, "frameshift_variant|splice_acceptor|splice_donnor|stop_gained")
         )
