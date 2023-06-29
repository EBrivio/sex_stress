# January 2020 - Elena Brivio
# scRNA-Seq on mouse PVN with 10x Genomics
# Combined analysis of two experiments:
# Controls vs Acute restraint (ARS) in M, F and OVX : ars.male.control or . stress
# Controls vs chronic mild stress (CMS) + ARS in M, F : cms.males.control or .stress

# Analysis on response distance

################### set-up environment   ###################

rm(list=ls())

library(Seurat)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(Matrix)
library(cowplot)
library(tidyverse)
library(AnnotationHub)
library(dplyr)
library(packcircles)
library(viridis)
library(ggrepel)
library(RColorBrewer)
library(RRHO)
library(RRHO2)
library(ComplexHeatmap)


sessionInfo()

#### Dataset import #####
pvn <- readRDS ("Objects/pvn_filtered_cleaned.Rds")

#### Differential gene expression import ####
# differential gene expression output
diff_gen <- read.csv(file = "differential_genes/combined.csv")
# sex = male or female
# exp = baseline or cms
# sex_exp = male_baseline, male_cms, female_baseline or female_cms
# cell_type = cell clusters of the dataset

#### N. DEGs import ####
# summary of DEGs per cell cluster
degs <- read.table("differential_genes/N_DEGs_bycluster_bysample_q0.05.txt", sep = "\t", header = T)

#### Clusters size import ####
# summary of cluster cell size
cluster_size <- read.table("Counts/counts_by_sample_cell.types.det.neu.txt", sep = "\t", header = T)

# Distance analysis CMS vs ARS ####
# Foldchange - n. genes 
fold_count_abs <- diff_gen %>% 
  dplyr::filter(p_val_BH < sig_value) %>% 
  dplyr::select(cell_type, gene, avg_logFC, sex_exp, exp) %>% 
  group_by(cell_type, sex_exp, exp) %>% 
  summarise(avg_logFC = median(abs(avg_logFC)),
            count = n()) 

fold_count_abs <- fold_count_abs %>% 
  left_join(cluster_size)

#logFC and count: geom distance
fold_count_abs_4plot <- fold_count_abs %>% 
  separate(sex_exp, into = c("sex", "exp"), remove = F) %>% 
  group_by(cell_type, sex) %>% 
  mutate(count1 = case_when(exp == "cms" ~ count, 
                            exp == "baseline" ~ -count),
         avg_logFC1 = case_when(exp == "cms" ~ abs(avg_logFC), 
                                exp == "baseline" ~ -abs(avg_logFC)),
         delta_N = sum(count1),
         delta_FC = sum(avg_logFC1),
         distance = sqrt((delta_N^2 + delta_FC^2))) %>% 
  group_by(sex) %>% 
  mutate(Z_distance = scale(distance))


# _Sex sum #####
fold_count_abs_4plot <- fold_count_abs_4plot %>% 
  group_by(cell_type) %>% 
  mutate(distance_MF = case_when(sex == "male" ~ Z_distance, 
                                 sex == "female" ~ 1/Z_distance),
         distance_FM = case_when(sex == "female" ~ Z_distance, 
                                 sex == "male" ~ 1/Z_distance),
         M_F_sum = sum(Z_distance)) %>%
  ungroup() %>% 
  group_by(measurement) %>% 
  mutate(Z_values = scale(values)) 