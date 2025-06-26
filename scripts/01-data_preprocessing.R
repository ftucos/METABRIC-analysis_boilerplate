library(genefu)
library(data.table)
library(tidyverse)

# Preprocess expression data  ------------------------------------------------------
exp_raw <- fread("data/raw_data/brca_metabric/data_mrna_illumina_microarray.txt") 
dup_symbol <- exp_raw %>% group_by(Hugo_Symbol) %>% summarize(n = n()) %>% filter(n > 1)

# expression data with unique symbols
exp_unique <- exp_raw %>%
  filter(!Hugo_Symbol %in% dup_symbol$Hugo_Symbol) %>%
  select(-Entrez_Gene_Id) %>%
  column_to_rownames("Hugo_Symbol") %>%
  as.matrix()

# expression data with duplicate symbols
# resolve by computing mean
exp_duplicate <- exp_raw %>% filter(Hugo_Symbol %in% dup_symbol$Hugo_Symbol) %>% select(-Entrez_Gene_Id) %>%
  group_by(Hugo_Symbol) %>%
  # alternative, select most variant probe
  summarize_all(~mean(., na.rm = T)) %>%
  column_to_rownames("Hugo_Symbol") %>%
  as.matrix()

# combine unique data
exp <- rbind(exp_unique, exp_duplicate)
# free some ram
rm(exp_unique, exp_duplicate)

# use genefu to obtain a complete THREEGENE/SCMOD2 classification for all the patients -----------
exp_entrez_raw <- fread("data/raw_data/brca_metabric/data_mrna_illumina_microarray.txt") %>%
  select(-Hugo_Symbol)

dup_ensembl_id <- exp_entrez_raw %>% group_by(Entrez_Gene_Id) %>% summarize(n = n()) %>% filter(n > 1)

# expression data with unique symbols
exp_entrez_unique <- exp_entrez_raw %>%
  filter(!Entrez_Gene_Id %in% dup_ensembl_id$Entrez_Gene_Id) %>%
  column_to_rownames("Entrez_Gene_Id") %>%
  as.matrix()

# expression data with duplicate symbols
# resolve by computing mean
exp_entrez_duplicate <- exp_entrez_raw %>%
  filter(Entrez_Gene_Id %in% dup_ensembl_id$Entrez_Gene_Id) %>%
  group_by(Entrez_Gene_Id) %>%
  # alternative, select most variant probe
  summarize_all(~mean(., na.rm = T)) %>%
  column_to_rownames("Entrez_Gene_Id") %>%
  as.matrix()

exp_entrez <- rbind(exp_entrez_unique, exp_entrez_duplicate)
# free some ram
rm(exp_entrez_unique, exp_entrez_duplicate)

annot <- data.frame(
  probe        = rownames(exp_entrez),    # required but can be any unique id
  #Gene.Symbol  = rownames(exp),    # <-- key column for mapping
  EntrezGene.ID  = rownames(exp_entrez),    # <-- key column for mapping
  stringsAsFactors = FALSE,
  row.names = rownames(exp_entrez)
)

data(scmod2.robust)

scmod2_classification <- molecular.subtyping(
  sbt.model = "scmod2",          # the pretrained model
  data      = t(exp_entrez),     # your expression matrix
  annot     = annot,         
  do.mapping = T
)

scmod2_class <- scmod2_classification[["subtype"]] %>% 
  as.data.frame() %>%
  rownames_to_column("PATIENT_ID") %>%
  rename("SCMOD2" =  ".")


# peprocess clinical ----------------------------
patient <- fread("data/raw_data/brca_metabric/data_clinical_patient.txt", skip = 4) %>%
  filter(PATIENT_ID %in% colnames(exp))

sample <- fread("data/raw_data/brca_metabric/data_clinical_sample.txt", skip = 4) %>%
  filter(PATIENT_ID %in% colnames(exp))

cBioportal_metadata <- full_join(patient, sample)

original_metadata <- rbind(
  fread("data/raw_data/brca_metabric_original_metadata/table_S2_revised.txt", na.strings = c("", "null")),
  fread("data/raw_data/brca_metabric_original_metadata/table_S3_revised.txt", na.strings = c("", "null"))
) %>%
  filter(METABRIC_ID %in% colnames(exp))

# 1980 patients --> removed 6 angiosarcoma and phylloides tumors --> 1974 patients
metadata <- full_join(
  original_metadata %>% rename(PATIENT_ID = METABRIC_ID),
  cBioportal_metadata %>% select(-NPI)) %>%
  left_join(scmod2_class) %>%
  select(PATIENT_ID, Sex = SEX, Age = AGE_AT_DIAGNOSIS, Iferred_Menopausal_State = INFERRED_MENOPAUSAL_STATE, Cohort = COHORT, NPI, Sample_tye = SAMPLE_TYPE,
         Cancer_type = CANCER_TYPE_DETAILED, Histotype = histological_type, Grade = GRADE, Tumor_size = size, Positive_nodes = LYMPH_NODES_EXAMINED_POSITIVE, Stage = TUMOR_STAGE,
         THREEGENE, SCMOD2, Pam50 = Pam50Subtype, Claudin_subtype = CLAUDIN_SUBTYPE, P53_mutation_status, TMB = TMB_NONSYNONYMOUS,
         Type_of_surgery = BREAST_SURGERY, Treatment, Chemotherapy = CHEMOTHERAPY, Hormonotherapy = HORMONE_THERAPY, Radiotherpay = RADIO_THERAPY, 
         RFS_MONTHS, RFS_STATUS, OS_MONTHS, OS_STATUS, VITAL_STATUS) %>%
  # remove Angiosarcoma and Phylloides tumors
  filter(Cancer_type != "Breast Angiosarcoma",
         !Histotype %in% c("PHYL"),
         # MB-0284 is indicated as benign but bears TP53 mutation
         ) %>%
  mutate(
    OS        = str_extract(OS_STATUS, "^[0-1]") %>% as.numeric(),
    OS_years  = OS_MONTHS/12,
    # Recode DSS from ex. vital status column
    DSS = case_when(
      VITAL_STATUS == "Living" ~ 0,
      VITAL_STATUS == "Died of Other Causes" ~ 0,
      VITAL_STATUS == "Died of Disease" ~ 1,
      TRUE ~ NA),
    DSS_years = OS_MONTHS/12,
    RFS       = str_extract(RFS_STATUS, "^[0-1]") %>% as.numeric(),
    RFS_years = RFS_MONTHS/12, .keep = "unused"
  ) %>%
  # note, not clear meaning of stage column, so we will not use it
  mutate(
    # more convenient scale to interpret HR in multivariable analysis
    Age_10_years_increase = Age/10,
    pT = case_when(
      Tumor_size == 0 ~ "Tis",
      # not using pT staging, but rather Tumor_size because no infomration for pT4 are available
      Tumor_size > 0 & Tumor_size <= 20 ~ "≤ 2 cm",
      Tumor_size > 2 & Tumor_size <= 5 ~ "> 2 - 5 cm",
      Tumor_size > 5 ~ "> 5 cm",
      TRUE ~ "TX") %>% factor(levels = c("≤ 2 cm",  "> 2 - 5 cm", "> 5 cm", "Tis", "TX")),
    pN = case_when(
      Positive_nodes == 0 ~ "N0",
      Positive_nodes > 0 & Positive_nodes <= 3  ~ "N1",
      Positive_nodes > 3 & Positive_nodes <= 9  ~ "N2",
      Positive_nodes > 9  ~ "N3",
      TRUE ~ "pNX") %>%
      factor(levels = c("N0", "N1", "N2", "N3"))) %>%
  relocate(pT, pN, .before = Tumor_size) %>%
  relocate(Age_10_years_increase, .after = Age)

# export ----------------
write_tsv(metadata, "data/processed/METABRIC-metadata.tsv")
saveRDS(metadata, "data/processed/METABRIC-metadata.rds")


exp_df <- exp %>%
  as.data.frame() %>%
  select(metadata$PATIENT_ID) %>%
  rownames_to_column("external_gene_name") 

exp_df_zscore <- exp_df %>%
  column_to_rownames("external_gene_name") %>%
  as.matrix() %>%
  t() %>% scale() %>% t() %>%
  as.data.frame() %>%
  rownames_to_column("external_gene_name") 

fwrite(exp_df, "data/processed/deduplicated_and_filtered-data_mrna_illumina_microarray.tsv.gz", compress = "gzip", sep = "\t", quote = FALSE)
fwrite(exp_df_zscore, "data/processed/deduplicated_and_filtered-data_mrna_illumina_microarray_zscores.tsv.gz", compress = "gzip", sep = "\t", quote = FALSE)

