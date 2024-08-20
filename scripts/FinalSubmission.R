# Sabin Hart
# August 29th
# Foundations of R - Final Submission


# libraries ---------------------------------------------------------------

library(tidyverse)
library(ggcorrplot)
library(ggthemes)
library(stringr)
library(rlist)


# data import -------------------------------------------------------------

setwd("/Users/sabinhart/Desktop/QBS 103 - R/final-proj")
gene_data_raw <- read_csv("./final-data/QBS103_GSE157103_genes.csv")
pheno_data_raw <- read_csv("./final-data/QBS103_GSE157103_series_matrix.csv")

gene_data <- gene_data_raw %>% 
  rename("gene_name" = 1) %>% 
  column_to_rownames("gene_name") %>% 
  t() %>% 
  data.frame() %>% 
  rownames_to_column("participant_id")

gene <- gene_data %>% 
  select(-1, -AADAC) %>% 
  select(sample(1:99, 15, replace = FALSE))

# combine data
cleaned_data <- gene_data %>% 
  merge(pheno_data_raw, by='participant_id') %>% 
  mutate(participant_id = as.numeric(sub(".*_(\\d+)_.*", "\\1", participant_id))) %>%
  mutate(participant_id = case_when(
    grepl('non', disease_status, ignore.case = TRUE) ~ participant_id + 200,
    TRUE ~ participant_id
  )) %>% 
  arrange(participant_id)

rm(gene_data, gene_data_raw, pheno_data_raw)

# correlation plot --------------------------------------------------------

corr <- round(cor(gene),2)

corr_plot <- ggcorrplot(corr, hc.order = TRUE, 
           type = "upper", 
           lab = TRUE, 
           lab_size = 2, 
           method="circle", 
           colors = c("firebrick1", "white", "green2"), 
           title="Gene Correlogram")

# ggsave("Correlogram.png", corr_plot, path = './plots')
rm(gene, corr, corr_plot)
