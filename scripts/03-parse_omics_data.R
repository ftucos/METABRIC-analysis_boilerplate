library(data.table)
library(tidyverse)
library(msigdbr)
library(glue)
library(ggstatsplot)
library(ggbeeswarm)
library(ggpubr)
library(ggh4x)

# Load Load Data
exp <- fread("data/processed/deduplicated_and_filtered-data_mrna_illumina_microarray.tsv.gz") %>%
  
metadata <- readRDS("data/processed/METABRIC-metadata.rds")

scores <- fread("data/interim/GSVA/singscore_scores.tsv")

# list genes of interest to be used for plots
exp_of_int <- c("CDK12")

gene_expression <- exp %>%
  filter(external_gene_name %in% exp_of_int) %>%
  column_to_rownames("external_gene_name") %>%
  t() %>%
  as.data.frame() %>%
  rename_all(~paste0(., "_EXP")) %>%
  rownames_to_column("PATIENT_ID")

CNA_of_int <- c("CDK12")

gene_CNA <- fread("data/brca_metabric/data_cna.txt") %>%
  select(-Entrez_Gene_Id) %>%
  filter(Hugo_Symbol %in% CNA_of_int) %>%
  column_to_rownames("Hugo_Symbol") %>%
  as.matrix() %>% t() %>% as.data.frame() %>%
  rename_all(~paste0(., "_CNA")) %>%
  mutate_all(as.character) %>%
  mutate_all(~recode(.,
    "-2" = "Deep/Homozygous Deletion",
    "-1" = "Shallow Deletion",
    "0" = "Diploid",
    "1" = "Gain",
    "2" = "Amplification"
  ) %>% factor(levels = c("Deep/Homozygous Deletion", "Shallow Deletion", "Diploid", "Gain", "Amplification"))) %>%
  rownames_to_column("PATIENT_ID")

metadata.1 <- metadata %>%
  full_join(gene_expression) %>%
  full_join(gene_CNA) %>%
  full_join(scores)

# used to plot all the scores at one with faceting variable
metadata.scores_long <-  metadata %>%
  full_join(gene_expression) %>%
  full_join(scores %>% gather(key="Signature", value = "GSVA", -1))

# plot --------------
ggplot(metadata.scores_long, aes(x=ABCC1_EXP, y=GSVA, group=Signature)) +
  geom_point(size = 0.3) +
  stat_smooth(method="lm", linewidth =0.7) +
  stat_cor(label.x = 6.5, label.y=-0.6, size = 8/.pt)+
  theme_bw(8) +
  facet_wrap(~Signature) +
  force_panelsizes(rows = unit(4.5, "cm"), cols = unit(5.5, "cm"))

ggsave("results/METABRIC/METABRIC-ABCC1_expression_and_stemness_signatures.pdf", width = 20, height = 20, unit = "cm" )

ggplot(metadata.scores_long %>% filter(Pam50 != "NC"), aes(x=ABCC1_EXP, y=GSVA, group=Signature)) +
  geom_point(size = 0.3) +
  stat_smooth(method="lm", linewidth =0.7) +
  stat_cor(label.x = 6.5, label.y=-0.6, size = 8/.pt)+
  theme_bw(8) +
  facet_grid(Pam50~Signature) +
  force_panelsizes(rows = unit(4.5, "cm"), cols = unit(5.5, "cm"))

ggsave("results/METABRIC/METABRIC-ABCC1_expression_and_stemness_signatures-by_subtype.pdf", width = 60, height = 30, unit = "cm" )


# ABCC1 expression and SRCIN1 CNA ---------
ggstatsplot::ggbetweenstats(metadata.1 %>% filter(Pam50 != "Her2"),
                            x=SRCIN1_CNA, y=ABCC1_EXP,
                            type = "parametric", bf.message = F, centrality.plotting = F)

metadata.1 %>% filter(Pam50 != "NC") %>%
  mutate(SRCIN1_CNA = factor(SRCIN1_CNA, levels = c("Diploid", "Deep/Homozygous Deletion", "Shallow Deletion", "Gain", "Amplification"))) %>%
  mutate(Pam50 = factor(Pam50, levels = c("LumA", "LumB", "Her2", "Basal", "Normal"))) %>%
  lm(ABCC1_EXP ~ SRCIN1_CNA + Pam50, data = .) %>%
    summary()
