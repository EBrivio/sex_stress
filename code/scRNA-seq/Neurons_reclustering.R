# Brivio et al.
# scRNA-Seq on mouse PVN with 10x Genomics
# Samples:
# Controls vs Acute restraint (ARS) in M, F and OVX : ars.male.control or . stress
# Controls vs chronic mild stress (CMS) + ARS in M, F : cms.males.control or .stress

# Neurons reclustering

rm(list = ls())

library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(fgsea)

sessionInfo()

# subset neurons ####
pvn <- readRDS ("Objects/pvn_filtered_cleaned.Rds")
neurons <- subset(pvn, cell.types.det == "Neurons")

# dataset set-up
neurons <- NormalizeData(object = neurons, 
                         normalization.method = "LogNormalize",
                         scale.factor = 10000)
neurons <- FindVariableFeatures(object = neurons, nfeatures = 4000, selection.method = 'vst', verbose = T)
length(x = VariableFeatures(object = neurons))

plot1 <- VariableFeaturePlot(neurons)
plot1 + coord_cartesian(ylim = c(0.5, 2)) 
top20 <- head(VariableFeatures(neurons), 20)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)

neurons <- ScaleData(object = neurons, vars.to.regress = NULL, verbose = T)
neurons <- RunPCA(object = neurons, pc.genes = neurons@var.genes, genes.print = 10)


#calculate Eigen values for PCs
eigValues_neu = (neurons$pca@stdev)^2  ## EigenValues

neurons <- JackStraw(object = neurons, dims = 40, num.replicate = 100)
neurons <- ScoreJackStraw(object = neurons, dims = 1:40)
ElbowPlot(neurons, ndims = 40)



# Approach: bootstrap analysis #####
# calculation of number of "stable clusters" and number of cells belonging to a stable cluster through a bootstrap analysis iterated across PCs, and resolution 

#PCs selection#####
##_PCs selection according to permutation ####
#random permutation (x50) of gene expression values across cell
data <- neurons@assays$RNA@data #lognorm version of counts
data <- as.matrix(data)
seurat <- neurons
eigValues <- data.frame(
                        eig = as.numeric(),
                        rep = as.numeric())

for (i in 1:50) {
  message(paste0("start of interation n.", i, " at ", Sys.time()))
  set.seed(i)
  data1 <- picante::randomizeMatrix(data, null.model = "richness")
  seurat@assays$RNA@data <- data1
  seurat <- ScaleData(object = seurat, vars.to.regress = NULL, verbose = F)
  seurat <- RunPCA(object = seurat, pc.genes = seurat@var.genes, verbose = F)
  eig = (seurat$pca@stdev)^2
  eig_df <- data.frame(eig = max(eig), rep = i)
  
  eigValues <- eigValues %>% bind_rows(eig_df)
  message(paste0("end of iteration n. ", i))
}

# saveRDS(eigValues, file = "PCs_test_eigValues.Rds")

print(eigValues)
ggplot(eigValues, aes(x = eig)) +
  geom_histogram(binwidth = 0.001) 
  # coord_cartesian(xlim = c(0,3))

mean(eigValues$eig)
pcs <- sum(eigValues_neu > mean(eigValues$eig)) #selection of PCs to use higher than mean eigenvalue of randomized dataset

## Test clusters stability#####

res <- seq(0.2, 2, 0.1)

source("function_stability_reclustering.R")

subsampling_clusters_stability(res = res, 
                               pcs = pcs, 
                               seurat_object = neurons,
                               iterations = 20, 
                               cells_subsampled = 0.5, 
                               output_dir = output_dir # define output directory
                                )

#Results bootstrap ####
pcs <- 38
results <- readRDS() # output from stability function


final_results_an <- results %>% 
  group_by(cluster, pcs, res, full_size) %>% 
  summarise(av_stability = mean(max)) %>% 
  mutate(stable = av_stability >= 0.5) %>%
  group_by(pcs, res)  %>% 
  mutate(tot_cl = max(cluster)+1)

#plot of stability as value 

#calculation number of stable clusters and number of cells in stable clusters
results_plot <- final_results_an %>% 
  filter(stable) %>% 
  group_by(pcs, res, tot_cl) %>% 
  mutate(n = n(),
         stable_cells = sum(full_size),
         perc = (n/tot_cl)*100) 

#export of these numbers into the parameters
#plots of the trend for each parameters: x= n. of cells in stable clusters, y = n. of stable clusters, fill = parameter that varies
ggplot(results_plot, aes(x = stable_cells, y = n, color = as.factor(res))) +
  geom_point() +
  labs(x = "Number of cells in stable clusters", y="Number of stable clusters")

  #identification of best paramaeter combination

final_results_an %>% 
  ggplot(aes(x = as.factor(res), y = av_stability)) +
  geom_point(aes(color = as.factor(cluster))) +
  geom_boxplot(alpha = 0.8) 

final_results_an %>% 
  select(-c(tot_cl)) %>% 
  distinct() %>% 
  ggplot(aes(x = as.factor(res), y = av_stability)) +
  # geom_bar() +
  stat_summary_bin(geom = "bar", fun = "mean") +
  geom_point(aes(color = as.factor(cluster)), alpha = 0.8) +
  stat_summary_bin(geom = "point", fun = "median", color = "red", pch = 3, size = 3) +
  geom_hline(yintercept = 0.5, color = "white", linetype = "dashed")
  

final_results_an %>% 
  group_by(res, pcs) %>% 
  summarise(stability = median(av_stability)) %>% 
  ggplot(aes(x = as.factor(res), y = stability)) +
  geom_point(alpha = 0.8) +
  ylab("Median Stability over clusters")


final_results_an %>% 
  group_by(res) %>% 
  summarise(m = median(av_stability)) %>% 
  ggplot(aes(x = as.factor(pcs), y = reorder(as.factor(res), m), fill = m)) +
  geom_tile() +
  theme_classic() +
  scale_fill_viridis_c() +
  ggtitle("median")


final_results_an %>% 
  group_by(res, pcs) %>% 
  summarise(m = median(av_stability)) %>% ungroup() %>%  top_n(n = 1, wt = m)

