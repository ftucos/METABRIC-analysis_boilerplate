library(data.table)
library(tidyverse)
library(msigdbr)
library(glue)
library(ggstatsplot)
library(ggbeeswarm)
library(ggpubr)
library(ggh4x)
library(survival)
library(survminer)
library(broom)
library(rlang)

# Load Load Data
exp <- fread("data/processed/deduplicated_and_filtered-data_mrna_illumina_microarray.tsv.gz")
  
metadata <- readRDS("data/processed/METABRIC-metadata.rds")

scores <- fread("data/interim/GSVA/singscore_scores.tsv")

# list genes of interest to be used for plots
exp_of_int <- c("MKI67", "ESR1", "PGR", "ERBB2")

gene_expression <- exp %>%
  filter(external_gene_name %in% exp_of_int) %>%
  column_to_rownames("external_gene_name") %>%
  t() %>%
  as.data.frame() %>%
  rename_all(~paste0(., "_EXP")) %>%
  rownames_to_column("PATIENT_ID")

CNA_of_int <- c("ERBB2")

gene_CNA <- fread("data/raw_data/brca_metabric/data_cna.txt") %>%
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
ggplot(metadata.1, aes(x=SCMOD2, y=MKI67_EXP, color=SCMOD2)) +
  geom_point(size = 0.3) +
  stat_smooth(method="lm", linewidth =0.7) +
  stat_cor(label.x = 6.5, label.y=0.1, size = 8/.pt)+
  theme_bw(8) +
  force_panelsizes(rows = unit(4.5, "cm"), cols = unit(5.5, "cm"))

ggsave("results/METABRIC/METABRIC-gene_expression_and_stemness_signatures.pdf", width = 20, height = 20, unit = "cm" )

ggplot(metadata.scores_long %>% filter(Pam50 != "NC"), aes(x=MKI67_EXP, y=GSVA, group=Signature)) +
  geom_point(size = 0.3) +
  stat_smooth(method="lm", linewidth =0.7) +
  stat_cor(label.x = 6.5, label.y=0.1, size = 8/.pt)+
  theme_bw(8) +
  facet_grid(Pam50~Signature) +
  force_panelsizes(rows = unit(4.5, "cm"), cols = unit(5.5, "cm"))

ggsave("results/METABRIC/METABRIC-gene_expression_and_stemness_signatures-by_subtype.pdf", width = 60, height = 30, unit = "cm" )

# ABCC1 expression and SRCIN1 CNA ---------
ggstatsplot::ggbetweenstats(metadata.1 %>% filter(Pam50 != "Her2"),
                            x=PGR_EXP, y=MKI67_EXP,
                            type = "parametric", bf.message = F, centrality.plotting = F)

metadata.1 %>% filter(Pam50 != "NC") %>%
  mutate(ERBB2_CNA = factor(ERBB2_CNA, levels = c("Diploid", "Deep/Homozygous Deletion", "Shallow Deletion", "Gain", "Amplification"))) %>%
  mutate(Pam50 = factor(Pam50, levels = c("LumA", "LumB", "Her2", "Basal", "Normal"))) %>%
  lm(ERBB2_EXP ~ ERBB2_CNA + Pam50, data = .) %>%
    summary()


# explore non-linear relationship with survival --------------------------------
cox_res <- coxph(Surv(OS_years, OS) ~ pspline(PGR_EXP), data = metadata.1)
cox_res %>% termplot(rug = TRUE, se = TRUE)

# explore multiple survival cutoff and thresholds ------------------------------

source("R/extract_cox_stats.R")
endpoint <- c("OS", "DSS", "RFS")
variable <- c("MKI67_EXP", "PGR_EXP")
quantiles <- list(0.1, 0.15, 0.25, 0.33, 0.5, 0.66, 0.75, 0.85, 0.9, "best")
# null is for all patients
subset <- list(NULL, 'SCMOD2 %in% c("ER+/HER2- Low Prolif", "ER+/HER2- High Prolif")', 'SCMOD2 == "HER2+"', 'SCMOD2 == "ER-/HER2-"')

# grid of model parameters to test
test_grid <- expand_grid(endpoint = endpoint,
                         variable = variable,
                         selected_quantile = quantiles,
                         .filter = subset)

# specify x as data for creating a partial function
map_function <- partial(extract_cox_stats, data = metadata.1)

res <- pmap_dfr(test_grid, map_function) %>%
  ungroup() %>%
  # highlight best quantile for each variable/endpoint combination
  mutate(best = ifelse(p.value == min(p.value), TRUE, FALSE), .by = c(variable, endpoint)) %>%
  # Manually adjust labels
  mutate(subgroup = ifelse(subgroup == "", "All", subgroup) %>%
           str_remove("^[^:]+: "))

# plot thresholds results
res %>%
  # add pvalue label
  mutate(signif = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01 ~ "**",
    p.value < 0.05 ~ "*",
    TRUE ~ "")) %>%
  ggplot(aes(y=as.character(quantile), x=endpoint, fill = HR)) +
  geom_tile() +
  geom_text(aes(label = signif)) +
  facet_wrap(~subgroup + variable, scales = "free") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1, na.value = "grey90") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  ylab("Quantile") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey90"))
