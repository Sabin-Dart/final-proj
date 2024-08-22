# Sabin Hart
# August 29th
# Foundations of R - Final Submission


# libraries ---------------------------------------------------------------

library(tidyverse)
library(ggcorrplot)
library(ggthemes)
library(stringr)
library(rlist)
library(stargazer)



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


# summary stats table -----------------------------------------------------

# 2 additional continuous (3 total) and 1 additional categorical variable (3 total)

# continuous: apacheii, ferritin(ng/ml), age

# categorical: stratify on sex; mechanical_ventilation, icu_status

categorical_table <- cleaned_data %>%
  select(sex, mechanical_ventilation, icu_status) %>% 
  filter(sex != 'unknown') %>% 
  group_by(sex) %>% 
  summarize(mv_n = sum(mechanical_ventilation=='yes'),
            mv_p = round(mean(mechanical_ventilation=='yes')*100, 1),
            icu_n = sum(icu_status == 'yes'),
            icu_p = round(mean(icu_status == 'yes')*100, 1),
            c = n()) %>% 
  mutate('Sex' = sex,
         `On Mechanical Ventilation` = paste0(mv_n, " (", mv_p, "%)"),
         `In ICU` = paste0(icu_n, " (", icu_p, "%)"),
         Count = c,
         .keep='none')

continuous_table <- cleaned_data %>%
  select(apacheii, `ferritin(ng/ml)`, age, sex) %>% 
  filter(sex != 'unknown') %>% 
  filter(apacheii != 'unknown') %>% 
  filter(`ferritin(ng/ml)` != 'unknown') %>% 
  mutate(age = as.numeric(age),
         apacheii = as.numeric(apacheii),
         `ferritin(ng/ml)` = as.numeric(`ferritin(ng/ml)`)) %>% 
  group_by(sex) %>% 
  summarize(a_mean = round(mean(age, na.rm = T), 1),
            a_sd = round(sd(age, na.rm = T), 1),
            ap_mean = round(mean(apacheii, na.rm = T), 1),
            ap_sd = round(sd(apacheii), 1),
            f_mean = round(mean(`ferritin(ng/ml)`, na.rm = T)),
            f_sd = round(sd(`ferritin(ng/ml)`))) %>% 
  mutate('Sex' = sex,
         Age = paste0(a_mean, " (", a_sd, ")"),
         `Apache II` = paste0(ap_mean, " (", ap_sd, ")"),
         `Ferritin (ng/mL)` = paste0(f_mean, " (", f_sd, ")"),
         .keep = 'none')

summary_table <- left_join(continuous_table, categorical_table, by = "Sex")
rm(categorical_table, continuous_table, table_data, test)
stargazer(summary_table,
          type = 'latex',
          title = "Summary Statistics",
          summary = FALSE,
          out = 'plots/summary_table.tex')


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
