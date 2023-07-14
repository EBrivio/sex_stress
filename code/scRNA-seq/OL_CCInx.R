# June 2020
# Elena Brivio
# Cell-Cell interaction analysis with CCInx package from https://github.com/BaderLab/CCInx
# Script adapted from Nagy et al. Nat. Neuro. 2020 (Supplementary_R_Script_2.R)


# Environment set-up ####
rm(list = ls())

# devtools::install_github("BaderLab/CCInx")
# devtools::install_github('Baderlab/scClustViz')
# BiocManager::install("DESeq")
# install.packages('varhandle')
# install.packages('abind')
library(CCInx)
library(scClustViz)
library(DESeq)
library(Seurat)
library(varhandle)
library(abind)
library(tidyr)
library(dplyr)
library(abind)
library(ggplot2)
library(ggrepel)
library(purrr)
library(RColorBrewer)
library(tibble)

sessionInfo()



output <- "Objects/Cell_Cell_interaction_rewritten/"
dir.create(file.path(getwd(), output))
output_results <- "Cell_cell_integration/rewritten/"
dir.create(file.path(getwd(), output_results))
output_images <- "Images/Filtered_nodoublets/Cell_cell_integration/rewritten/final/"
dir.create(file.path(getwd(), output_images), recursive = T)

# Import data ####
pvn <- readRDS("Objects/pvn_filtered_cleaned.Rds")

# remove "_" in cell names because No "_" allowed in CCInx
pvn$cell.types.det.neu_ <- plyr::mapvalues(x = pvn$cell.types.det.neu, from = c("Neurons_mixed", "Mixed_Endothelial"), to = c("NeuronsMixed", "MixedEndothelial"))

cells1 <- "Oligodendrocytes"
cells2 <- c("Glut")
# cells2 <- c("Glut", "GABA", "NeuronsMixed")

# Create a for loop here if multiple cells2 want to be tested

# Subset dataset ####
o_g <- subset(pvn, cell.types.det.neu_ == cells1 | cell.types.det.neu_ == cells2)
table(o_g$cell.types.det.neu_, o_g$orig.ident)
table(o_g$cell.types.det.neu_)/8

pvn_test <- o_g

tt <- table(pvn_test$cell.types.det.neu_)/2
tt

pvn_test@meta.data$cells <- rownames(pvn_test@meta.data)
pvn_test@meta.data <- pvn_test@meta.data %>% 
  separate(cells, sep = "_", into = c("trash", "cells")) %>% 
  select(-trash)
pvn_test@meta.data$condition <- ""
# rm(pvn, o_g)

# Filter genes #####
# Determine the expression threshold based on markers expression levels 
FeaturePlot(pvn_test, "Mag", sort.cell = T)
FeaturePlot(pvn_test, c("Mag", "Gad1"), sort.cell = T)
FeaturePlot(pvn_test, c("Mag", "Gad1", "Plp1"), sort.cell = T)
FeaturePlot(pvn_test, c("Mag", "Gad1", "Plp1"), sort.cell = T, min.cutoff = 2.75)
FeaturePlot(pvn_test, c("Gad1"), sort.cell = T, min.cutoff = 2.75) #not ok for neurons! too stringent
FeaturePlot(pvn_test, c("Slc17a6"), sort.cell = T, min.cutoff = 1.5) 
FeaturePlot(pvn_test, c("Map2"), sort.cell = T, min.cutoff = 1.5) 
FeaturePlot(pvn_test, c("Avp"), sort.cell = T, min.cutoff = 1.5) 

FeaturePlot(pvn_test, c("Gad1"), sort.cell = T, min.cutoff = 2)
FeaturePlot(pvn_test, c("Gad1"), sort.cell = T, min.cutoff = 1.5)
FeaturePlot(pvn_test, c("Gad1"), sort.cell = T, max.cutoff = 1.5)
FeaturePlot(pvn_test, c("Mag"), sort.cell = T, min.cutoff = 2.75) #ok for oligos
FeaturePlot(pvn_test, c("Gad2"), sort.cell = T) 
FeaturePlot(pvn_test, c("Gad2"), sort.cell = T, min.cutoff = 1.5) 
FeaturePlot(pvn_test, c("Mag", "Gad1", "Plp1"), sort.cell = T, min.cutoff = 1.5, slot = "scale.data")


# Permutation Analysis####
#loop n_iter times
n_iter = 100
threshold_expression <- 1.5

two <- read.csv(text = "pair, Gene_A, cell_A, Gene_B, cell_B, condition, edgeWeight, diff_edges, repetition, exp_sex")


for (i in 1:n_iter){
  #subsetted into groups of cells mixed from both samples -> the new groups are neither controls nor stress
   ind <- c(sample(unique(pvn_test@meta.data[pvn_test$cell.types.det.neu_ == cells2, "cells"]), tt[cells2]), sample(unique(pvn_test@meta.data[pvn_test$cell.types.det.neu_ == cells1, "cells"]), tt[cells1])) #1/2 of the dataset
  #GABA 1327
  #Glut 713
  #Neuronsmixed
  #AVP 471.5
  
  groupA_cells <- which(pvn_test@meta.data[, "cells"] %in% ind) 
  groupB_cells <- which(!(pvn_test@meta.data[, "cells"] %in% ind)) #I select 1/2 of the dataset vs 1/2 of the dataset
  
  #replace condition meta data with new group labels
  pvn_test@meta.data[c(groupA_cells),]$condition <- c("groupA")
  pvn_test@meta.data[c(groupB_cells),]$condition <- c("groupB")
  
  #create cluster*condition meta data labels
pvn_test@meta.data$cluster_condition <- paste(pvn_test$cell.types.det.neu_, pvn_test$condition, sep = ".")


gsl <- BuildGeneStatList(inD = pvn_test,
                         cl = pvn_test@meta.data$cluster_condition,
                         assayType = "RNA",
                         assaySlot = "data") #"logcounts
inx <- BuildCCInx(GeneStatList = gsl,
                  Species = "mmusculus")


#run functions from CCinx filtering the connections
  connectA <- FilterInx_step1(inx, paste(cells2, 'groupA', sep = "."), paste(cells1, 'groupA', sep = "."), 'Ligand', 'Receptor')
  connectA <- FilterInx_GeneMagnitude(connectA, threshold_expression)
  connectA <- FilterInx_genenames(connectA, connectA$nodes$gene)
  
  connectB <- FilterInx_step1(inx, paste(cells2, 'groupB', sep = "."), paste(cells1, 'groupB', sep = "."), 'Ligand', 'Receptor')
  connectB <- FilterInx_GeneMagnitude(connectB, threshold_expression)
  connectB <- FilterInx_genenames(connectB, connectB$nodes$gene)
  
  message("Got lists direction 1; now merging 1...\n")
  
  AA <- connectA$edges
  # head(AA)
  AA <- AA %>% 
    separate(nodeA, into = c("Gene_A", "cell_A", "node"), sep = "_") %>% 
    separate(cell_A, into = c("cell_A", "condition")) %>% 
    select(-node) %>% 
    separate(nodeB, into = c("Gene_B", "cell_B", "node"), sep = "_") %>% 
    separate(cell_B, into = c("cell_B", "condition")) %>% 
    select(-node) %>% 
    mutate(pair = paste(Gene_A, cell_A, Gene_B, cell_B, sep = "_"))
  head(AA)
  
  BB <- connectB$edges
  head(BB)
  BB <- BB %>% 
    separate(nodeA, into = c("Gene_A", "cell_A", "node"), sep = "_") %>% 
    separate(cell_A, into = c("cell_A", "condition")) %>% 
    select(-node) %>% 
    separate(nodeB, into = c("Gene_B", "cell_B", "node"), sep = "_") %>% 
    separate(cell_B, into = c("cell_B", "condition")) %>% 
    select(-node) %>% 
    mutate(pair = paste(Gene_A, cell_A, Gene_B, cell_B, sep = "_"))
  
  AB <- AA %>%
    bind_rows(BB) %>% 
    group_by(pair) %>% 
    nest() %>% 
    mutate(
      diff_edges = purrr::map_dbl(data, function(d){
        if (length(d$edgeWeight) == 2) {
          p <- pivot_wider(data = d, id_cols = c(condition, edgeWeight),
                           names_from = condition, values_from = edgeWeight)
          p$groupB - p$groupA} else {
            c(-d[d$condition == "groupA",]$edgeWeight,
              d[d$condition == "groupB",]$edgeWeight)}
      })
    ) %>% 
    unnest(cols = data) %>% 
    mutate(exp_sex = "permutation",
           repetition = i)
  
  
  # ggplot(AA, aes(x = edgeWeight)) +
  #   geom_histogram()
  # ggplot(BB, aes(x = edgeWeight)) +
  #   geom_histogram()
  # ggplot(inx$edges, aes(x = edgeWeight)) +
  #   geom_histogram()
  # plots QC ###
  # cowplot::plot_grid(ggplot(AA, aes(x = edgeWeight)) +
  #                      geom_histogram(),
  #                    ggplot(BB, aes(x = edgeWeight)) +
  #                      geom_histogram(),
  #                    ncol = 1)
  # ggplot(AA, aes(x = edgeWeight)) +
  #   geom_histogram() +
  #   xlim(0,1)
  # ggplot(BB, aes(x = edgeWeight)) +
  #   geom_histogram() +
  #   xlim(0,1)
  # 
  # cowplot::plot_grid(ggplot(AA, aes(x = edgeWeight)) +
  #                      geom_histogram()+
  #                      xlim(0,1),
  #                    ggplot(BB, aes(x = edgeWeight)) +
  #                      geom_histogram()+
  #                      xlim(0,1),
  #                    ncol = 1)
  # # 
  # ggplot(AB, aes(x = diff_edges)) +
  #   geom_histogram(bins = 50) +
  #   xlim(-0.2, 0.2)
  # 
  # head(AB)
  # ggplot(AB, aes(x = pair, y = edgeWeight, fill = condition)) +
  #   geom_point(pch = 21) +
  #   coord_flip() +
  #   theme_classic() #looks ok

  AB %>%
    ungroup() %>%
    sample_n(size = 50, replace = T) %>%
    ggplot(aes(x = pair, y = edgeWeight, fill = condition)) +
      geom_point(pch = 21, size = 4) +
      coord_flip() +
      theme_classic() +
      theme(panel.grid.major.y = element_line(color = "grey", size = 0.5))

  message("Done direction 1.\n")
  
  message("Now starting direction 2...\n")
  
  connectA <- FilterInx_step1(inx, paste(cells1, 'groupA', sep = "."), paste(cells2, 'groupA', sep = "."), 'Ligand', 'Receptor')
  connectA <- FilterInx_GeneMagnitude(connectA, threshold_expression)
  connectA <- FilterInx_genenames(connectA, connectA$nodes$gene)
  
  connectB <- FilterInx_step1(inx, paste(cells1, 'groupB', sep = "."), paste(cells2, 'groupB', sep = "."), 'Ligand', 'Receptor')
  connectB <- FilterInx_GeneMagnitude(connectB, threshold_expression)
  connectB <- FilterInx_genenames(connectB, connectB$nodes$gene)
  
  message("Got lists direction 2; now merging 2...\n")
  
  AA <- connectA$edges
  # head(AA)
  AA <- AA %>% 
    separate(nodeA, into = c("Gene_A", "cell_A", "node"), sep = "_") %>% 
    separate(cell_A, into = c("cell_A", "condition")) %>% 
    select(-node) %>% 
    separate(nodeB, into = c("Gene_B", "cell_B", "node"), sep = "_") %>% 
    separate(cell_B, into = c("cell_B", "condition")) %>% 
    select(-node) %>% 
    mutate(pair = paste(Gene_A, cell_A, Gene_B, cell_B, sep = "_"))
  head(AA)
  
  BB <- connectB$edges
  head(BB)
  BB <- BB %>% 
    separate(nodeA, into = c("Gene_A", "cell_A", "node"), sep = "_") %>% 
    separate(cell_A, into = c("cell_A", "condition")) %>% 
    select(-node) %>% 
    separate(nodeB, into = c("Gene_B", "cell_B", "node"), sep = "_") %>% 
    separate(cell_B, into = c("cell_B", "condition")) %>% 
    select(-node) %>% 
    mutate(pair = paste(Gene_A, cell_A, Gene_B, cell_B, sep = "_"))
  
  BA <- AA %>%
    bind_rows(BB) %>% 
    group_by(pair) %>% 
    nest() %>% 
    mutate(
      diff_edges = purrr::map_dbl(data, function(d){
        if (length(d$edgeWeight) == 2) {
          p <- pivot_wider(data = d, id_cols = c(condition, edgeWeight), names_from = condition, values_from = edgeWeight)
          p$groupB - p$groupA} else {
            c(-d[d$condition == "groupA",]$edgeWeight,
              d[d$condition == "groupB",]$edgeWeight)}
      })
    ) %>% 
    unnest(cols = data) %>% 
    mutate(exp_sex = "permutation",
           repetition = i)
  head(BA)
  two <- two %>%
    bind_rows(AB, BA)
  message("Done Direction 2.\n")
  message(paste0(i, "th iteration.\n"))
  
  if (i == n_iter){
    print(
      two %>% 
        mutate(cell_AB = paste(cell_A, cell_B, sep = ".")) %>% 
        group_by(repetition, cell_AB) %>% 
        summarise(n = n()))
  }
}

save.image(file = paste0(output, "Cell_cell_interaction_", cells1, "_", cells2, "_", threshold_expression, "_iterations.RData"))

# head(two)
# tail(two)
write.csv(two, file = paste0(output_results, "Cell_cell_interaction_", cells1, "_", cells2, "_", threshold_expression, "_iterations.csv"), row.names = F)

two <- read.csv(file = paste0(output_results, "Cell_cell_interaction_", cells1, "_", cells2, "_", threshold_expression, "_iterations.csv"))

