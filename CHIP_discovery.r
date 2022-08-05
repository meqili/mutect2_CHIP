library(openxlsx)
library(data.table)
library(tidyverse)
library(stringi)

# if two files was not read in environment
CHIP_mutations <- read.xlsx("D:/2022-04_Mutect/mutect_CHIP/Supplementary Table 2. Leukemogenic driver mutations queried in this study.xlsx")
CHIP_gene_coordinates <- fread("D:/2022-04_Mutect/mutect_CHIP/CHIP_gene_coordinates_ENST.tsv")

CHIP_mutations <- left_join(CHIP_mutations, CHIP_gene_coordinates, 
                            by = c("Accession" = "refseq_mrna"))

#testing <- paste(CHIP_mutations$`Reported.mutations.used.for.variant.calling.-.ql2387`, collapse = ", ")
#testing <- str_split(testing, ",\\s*" )[[1]]

CHIP.variants.list <- c()

for (ENST_ID in CHIP_mutations$ensembl_transcript_id) {
  chip_mutations <- CHIP_mutations[CHIP_mutations$ensembl_transcript_id == ENST_ID,]$`Reported.mutations.used.for.variant.calling.-.ql2387`
  chip_mutations <- str_split(chip_mutations, ",\\s*" )[[1]]

  ## 3rd type
  ## those two needed to check manually
  ## FY590-591GD (FLT3), W515-518KT (MPL)
  AA_change_list <- chip_mutations[which(str_detect(chip_mutations, "[:upper:]\\d+") &
                                           !str_detect(chip_mutations, "-"))]
  if (length(AA_change_list) >0) {
    CHIP.variants.list <- variants.list %>% filter(Feature_ID == ENST_ID) %>%
      filter(new.HGVS.p %in% AA_change_list) %>%
      rbind(CHIP.variants.list)}
  
  ## 2nd type
  Deletion_insertion_format <- regex("
  (del/ins|del|ins)     # replacement type
  (\\d+\\-\\d+|\\d+)         # replace posistion
  ([A-Z]+)?             # replaced protein
  ", comments = TRUE)
  
  delins.p <- chip_mutations[which(str_detect(chip_mutations, "^del|^ins"))]
  if (length(delins.p) >0) {
    delins.p.pattern <- str_match(delins.p, Deletion_insertion_format)
    delins.p.format <- paste0(gsub("NA", "", paste0("[A-Z]", 
                                                    gsub("-", "_[A-Z]", delins.p.pattern[,3]),
                                                    gsub("/", "", delins.p.pattern[,2]),
                                                    delins.p.pattern[,4])), collapse = "|")
    CHIP.variants.list <- variants.list %>% filter(Feature_ID == ENST_ID) %>%
      filter(str_detect(new.HGVS.p, delins.p.format)) %>%
      rbind(CHIP.variants.list)}
  
  ## 1st type
  general_types <- chip_mutations[
    which(str_detect(chip_mutations, "Frameshift|nonsense|splice-site|missense"))]
  if (length(general_types) >0) {
    for (general_type in general_types) {
      general_pattern <- str_replace_all(
        str_split(general_type, "\\s+" )[[1]][1], "/", "|")
      general_pattern <- stri_replace_all_fixed(general_pattern,
                                                pattern = c("Frameshift", "nonsense", "splice-site", "missense"),
                                                replacement = c("frameshift_variant", "stop_gained", "splice_acceptor|splice_donnor", "missense_variant"),
                                                vectorize_all = FALSE)
      
      general_loc <- str_split(general_type, "\\s+" )[[1]][3]
      if (! is.na(general_loc) & general_loc == "protein") {
        column_to_filter <- "AA.pos/AA.length"
      } else if (! is.na(general_loc) & general_loc == "cDNA") {
        column_to_filter <- "cDNA.pos/cDNA.length"
      }else  { column_to_filter <- "Rank" }
      
      general_pos <- str_split(general_type, "\\s+" )[[1]][4]
      if (is.na(general_pos)) {
        general_pos_pattern <- "."
      } else {
        general_pos_pattern <- paste0(paste0("^",
                                             c(as.numeric(str_split(general_pos, "-")[[1]])[1]:as.numeric(str_split(general_pos, "-")[[1]])[2]), 
                                             "/"), collapse = "|")
      }
      
      CHIP.variants.list <- variants.list %>% filter(Feature_ID == ENST_ID) %>%
        filter(str_detect(get(column_to_filter), general_pos_pattern) & 
                 str_detect(Annotation, general_pattern)) %>%
        rbind(CHIP.variants.list)
    }}}

## FY590-591GD (FLT3), W515-518KT (MPL)
ENST_ID <- CHIP_mutations[CHIP_mutations$Gene.name == "FLT3",]$ensembl_transcript_id
general_pos_pattern <- paste0(paste0("^", c(590,591), 
                                     "/"), collapse = "|")
variants.list %>% filter(Feature_ID == ENST_ID) %>%
  filter(str_detect(`AA.pos/AA.length`, general_pos_pattern)) %>% 
  write.xlsx("FLT3_FY590-591GD_regions.xlsx")

ENST_ID <- CHIP_mutations[CHIP_mutations$Gene.name == "MPL",]$ensembl_transcript_id
general_pos_pattern <- paste0(paste0("^", c(512:520), 
                                     "/"), collapse = "|")
variants.list %>% filter(Feature_ID == ENST_ID) %>%
  filter(str_detect(`AA.pos/AA.length`, general_pos_pattern)) %>% 
  write.xlsx("MPL_W515-518KT_regions.xlsx")

## add Telo#########
all_caseData_qPCR <- read.xlsx("D:/2022-04_Mutect/2022-05-31_final_cleaned_telomere_length_ql2387.xlsx") %>%
  select(sample_internal_name, percentile)
# CHIP.variants.list <- read.xlsx("CHIP_variants_foundin_allIPFsamples.xlsx")
CHIP.variants.list <- left_join(CHIP.variants.list, all_caseData_qPCR, 
                                by=c("Indiv" = "sample_internal_name"))
###### final #########
write.xlsx(CHIP.variants.list, file = "CHIP_variants_foundin_all_IPFsamples_telo.xlsx")
###################


### visu
CHIP.variants.list %>% select(Indiv, Gene_Name) %>% 
  unique() %>% select(Gene_Name) %>% 
  table() %>% sort() %>% 
  barplot(las=2, cex.names=.7, xlab="Genes",
          ylab="Number of samples have at least one CHIP mutations")
  
CHIP.variants.list %>% filter(Gene_Name == "STAG2") %>% 
  ggplot(aes(x = POS)) +
  geom_histogram(binwidth = 1)

CHIP.variants.list %>% filter(Gene_Name == "STAG2" & POS >=123184949  & POS<= 123184989  ) %>%
  ggplot(aes(x = POS)) + geom_histogram(binwidth = 1) +xlab("") + theme(axis.text.x=element_blank())
