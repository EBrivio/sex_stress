# Brivio et al.
# scRNA-Seq on mouse PVN with 10x Genomics
# Samples:
# Controls vs Acute restraint (ARS) in M, F and OVX : ars.male.control or . stress
# Controls vs chronic mild stress (CMS) + ARS in M, F : cms.males.control or .stress

# Analysis of ambient RNA

#### Set-up environment ####
rm(list = ls())

library(Seurat)
library(ggplot2)
library(tidyverse)

sessionInfo()


### Importing empty droplets ####
# minimum number of cells containing a gene in order to be kept
minimum.cells = 3

### Loading data from the cell ranger output files + generate a Seurat object per each sample
# ARS
gg <- Read10X(data.dir = "Raw_matrices/Acute/count_data/mpg_L15893_F-c_S13/raw_feature_bc_matrix/")
ars.female.control <- CreateSeuratObject(counts = gg, project = "ars.female.control", min.cells = minimum.cells)
afc <- FeatureScatter(object = ars.female.control, feature2 = "nFeature_RNA", feature1 = "nCount_RNA", pt.size = 0) + NoLegend()

gg <- Read10X(data.dir = "Raw_matrices/Acute/count_data/mpg_L15893_F-s_S14/raw_feature_bc_matrix/")
ars.female.stress <- CreateSeuratObject(counts = gg, project = "ars.female.stress", min.cells = minimum.cells)
afs <-FeatureScatter(object = ars.female.stress, feature2 = "nFeature_RNA", feature1 = "nCount_RNA", pt.size = 0) + NoLegend()

gg <- Read10X(data.dir = "Raw_matrices/Acute/count_data/mpg_L15893_M-c_S15/raw_feature_bc_matrix/")
ars.male.control <- CreateSeuratObject(counts = gg, project = "ars.male.control", min.cells = minimum.cells)
amc <-FeatureScatter(object = ars.male.control, feature2 = "nFeature_RNA", feature1 = "nCount_RNA", pt.size = 0) + NoLegend()

gg <- Read10X(data.dir = "Raw_matrices/Acute/count_data/mpg_L15893_M-s_S16/raw_feature_bc_matrix/")
ars.male.stress <- CreateSeuratObject(counts = gg, project = "ars.male.stress", min.cells = minimum.cells)
ams <-FeatureScatter(object = ars.male.stress, feature2 = "nFeature_RNA", feature1 = "nCount_RNA", pt.size = 0) + NoLegend()


# CMS + ARS
gg <- Read10X(data.dir = "Raw_matrices/Chronic/count_matrices/female_control/raw_feature_bc_matrix/")
cms.female.control <- CreateSeuratObject(counts = gg, project = "cms.female.control", min.cells = minimum.cells)
cfc <- FeatureScatter(object = cms.female.control, feature2 = "nFeature_RNA", feature1 = "nCount_RNA", pt.size = 0) + NoLegend()


gg <- Read10X(data.dir = "Raw_matrices/Chronic/count_matrices/female_stress/raw_feature_bc_matrix/")
cms.female.stress <- CreateSeuratObject(counts = gg, project = "cms.female.stress", min.cells = minimum.cells)
cfs <- FeatureScatter(object = cms.female.stress, feature2 = "nFeature_RNA", feature1 = "nCount_RNA", pt.size = 0) + NoLegend()

gg <- Read10X(data.dir = "Raw_matrices/Chronic/count_matrices/male_control/raw_feature_bc_matrix/")
cms.male.control <- CreateSeuratObject(counts = gg, project = "cms.male.control", min.cells = minimum.cells)
cmc <- FeatureScatter(object = cms.male.control, feature2 = "nFeature_RNA", feature1 = "nCount_RNA", pt.size = 0) + NoLegend()

gg <- Read10X(data.dir = "Raw_matrices/Chronic/count_matrices/male_stress/raw_feature_bc_matrix/")
cms.male.stress <- CreateSeuratObject(counts = gg, project = "cms.male.stress", min.cells = minimum.cells)
cms <- FeatureScatter(object = cms.male.stress, feature2 = "nFeature_RNA", feature1 = "nCount_RNA", pt.size = 0) + NoLegend()

rm(gg)

cowplot::plot_grid(afc, afs, amc, ams, cfc, cfs, cmc, cms, nrow = 2)

# Filtering single objects ####

## Create one Seurat object combining all samples together ####
empty <- merge(x = ars.female.control, y = c(ars.female.stress, ars.male.control, ars.male.stress, cms.female.control, cms.female.stress, cms.male.control, cms.male.stress), add.cell.ids = c("ars.female.control", "ars.female.stress", "ars.male.control", "ars.male.stress", "cms.female.control", "cms.female.stress", "cms.male.control", "cms.male.stress"))

# Features
ggplot(empty@meta.data, aes(empty@meta.data$nFeature_RNA)) +
  geom_bar() +
  theme_minimal() + 
  coord_cartesian(ylim = c(0,1000), xlim = c(0,1000)) 
  # scale_x_continuous(breaks = seq(0, 300, 20), 
  #                    minor_breaks = seq(0, 300, 5))

ggplot(empty@meta.data, aes(empty@meta.data$nFeature_RNA)) +
  geom_bar() +
  theme_minimal() + 
  coord_cartesian(ylim = c(0,2000), xlim= c(35, 140)) +
  scale_x_continuous(breaks = seq(0, 300, 20), 
                     minor_breaks = seq(0, 300, 5))

empty = subset(x = empty, subset = nFeature_RNA < 140 & nFeature_RNA >35) 

# UMI
ggplot(empty@meta.data, aes(empty@meta.data$nCount_RNA)) +
  geom_bar() +
  theme_minimal() + 
  coord_cartesian(ylim = c(0,2000)) +
  scale_x_continuous(breaks = seq(0, 300, 20), minor_breaks = seq(0, 300, 5))

ggplot(empty@meta.data, aes(empty@meta.data$nCount_RNA)) +
  geom_bar() +
  theme_minimal() + 
  coord_cartesian(ylim = c(0,2000), xlim= c(0, 150)) +
  scale_x_continuous(breaks = seq(0, 300, 20), minor_breaks = seq(0, 300, 5))

empty = subset(x = empty, subset = nCount_RNA < 150) 

ggplot(empty@meta.data, aes(empty@meta.data$nCount_RNA)) +
  geom_bar() +
  theme_minimal()
ggplot(empty@meta.data, aes(empty@meta.data$nFeature_RNA)) +
  geom_bar() +
  theme_minimal()



VlnPlot(object = empty, features = c("nFeature_RNA", "nCount_RNA"), group.by="orig.ident", pt.size = 0) + NoLegend()

FeatureScatter(object = empty, feature1 = "nCount_RNA", feature2  = "nFeature_RNA", group.by = "orig.ident")


rm(ars.female.control, ars.female.stress, ars.male.control, ars.male.stress, cms.male.control, cms.male.stress, cms.female.control, cms.female.stress)

# Log normalization

empty <- NormalizeData(object = empty, normalization.method = "LogNormalize", 
                     scale.factor = 10000)


empty <- FindVariableFeatures(object = empty, nfeatures = 4000, selection.method = 'vst', verbose = T)
VariableFeaturePlot(empty) + coord_cartesian(ylim = c(0.99, 1.1)) 
empty <- ScaleData(object = empty, vars.to.regress = NULL, verbose = T)

# Calculation Ambient ####
FeatureScatter(empty, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident" )

# Average expression ambient
empty$all <- "none"
Idents(empty) <- "all"
avg <- log1p(AverageExpression(empty, verbose = T, slot = "data")$RNA)
avg$gene <- rownames(avg)
avg <- avg %>% rename("expression" = "none")


ggplot(avg, aes(x = expression)) +
  geom_histogram(binwidth = 0.05) +
  theme_minimal()

# top 50
avg %>% 
  arrange(-expression) %>% 
  top_n(n = 50, wt = expression) %>% 
ggplot(aes(y = reorder(gene, expression), x = expression)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_x_continuous(expand = c(0,0), breaks = seq(0, 5, 0.5)) +
  geom_text(aes(label = gene, x = 0, y = gene), hjust = 0) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())

avg <- avg %>% mutate(droplet = "empty")


# Bar plots of GOI
pvn <- readRDS(file = "Objects/pvn_filtered_cleaned.Rds") # full cleaned dataset
goi <- c("Avp", "Oxt", "Crh", "Sst", "Trh")
avg %>% 
  filter(gene %in%goi) %>% 
  ggplot(aes(x = gene, y = expression)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_y_continuous(expand = c(0,0), breaks = seq(0, 5, 0.5)) +
  ggtitle("Ambient RNA")+
  theme(text = element_text(family = "Helvetica", size = 20), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5)) 

