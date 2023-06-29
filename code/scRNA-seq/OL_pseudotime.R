# Brivio et al.
# scRNA-Seq on mouse PVN with 10x Genomics
# Samples:
# Controls vs Acute restraint (ARS) in M, F and OVX : ars.male.control or . stress
# Controls vs chronic mild stress (CMS) + ARS in M, F : cms.males.control or .stress

# Analysis of pseudotime for oligodendrocytes

### Set-up ####
rm(list = ls())
library(dplyr)
library(tidyr)
library(Seurat)
library(ggplot2)
library(monocle3)
library(report)


sessionInfo()

# Subset oligodendrocytes ####
pvn <- readRDS ("Objects/pvn_filtered_cleaned.Rds")
oligos <- subset(pvn, cell.types == "Oligodendrocytes")

## Reclustering oligodendrocytes ####

oligos <- NormalizeData(object = oligos, normalization.method = "LogNormalize", scale.factor = 10000)
oligos <- FindVariableFeatures(object = oligos, nfeatures = 4000, selection.method = 'vst', verbose = T)
length(x = VariableFeatures(object = oligos))

plot1 <- VariableFeaturePlot(oligos) #inspect variable genes
plot1 + coord_cartesian(ylim = c(0.5, 2))
top20 <- head(VariableFeatures(oligos), 20)
plot2 <- LabelPoints(plot = oligos, points = top20, repel = TRUE)

oligos <- ScaleData(object = oligos, vars.to.regress = NULL, verbose = T)

oligos <- RunPCA(object = oligos, pc.genes = oligos@var.genes, genes.print = 10)

# _Clustering: 15PC ####
pcs <- 15
res <- 0.6
oligos <- FindNeighbors(object = oligos, dims = 1:pcs)
oligos <- FindClusters(object = oligos, reduction.type = "pca", resolution = res, verbose = T, save.SNN = TRUE, force.recalc =TRUE)
oligos <- RunUMAP(oligos, dims = 1:pcs)


# Creation of Monocle3 object ####

gene_annotation <- as.data.frame(rownames(oligos@reductions[["pca"]]@feature.loadings), row.names = rownames(oligos@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"

# part two, cell information
cell_metadata <- as.data.frame(oligos@assays[["RNA"]]@counts@Dimnames[[2]], row.names = oligos@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"

# part three, counts sparse matrix
New_matrix <- oligos@assays[["RNA"]]@counts
New_matrix <- New_matrix[rownames(oligos@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix

### Construct the basic cds object
oligos_mono <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)
# _transfer of information from Seurat to Monocle 3####
### Construct and assign the made up partition
recreate.partition <- c(rep(1, length(oligos_mono@colData@rownames)))
names(recreate.partition) <- oligos_mono@colData@rownames
recreate.partition <- as.factor(recreate.partition)

oligos_mono@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

### Assign the cluster info
list_cluster <- oligos@meta.data[["named.clusters"]]
names(list_cluster) <- oligos@assays[["RNA"]]@data@Dimnames[[2]]

oligos_mono@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
oligos_mono@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"

### Assign UMAP coordinate
oligos_mono@int_colData@listData[["reducedDims"]][["UMAP"]]<- oligos@reductions[["umap"]]@cell.embeddings

### Assign feature loading for downstream module analysis
oligos_mono@preprocess_aux$gene_loadings <- oligos@reductions[["pca"]]@feature.loadings


# Trajectory analysis ####
# partitioning
oligos_mono <- learn_graph(oligos_mono, use_partition = T)
plot_cells(oligos_mono, color_cells_by = "partition")

oligos_mono@colData$orig.ident <- oligos@meta.data$orig.ident

#trajectory graph
plot_cells(oligos_mono,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

# identify root on the plot 
rownames(oligos_mono@principal_graph_aux[['UMAP']]$dp_mst) <- NULL
colnames(oligos_mono@int_colData@listData[["reducedDims"]][["UMAP"]]) <- NULL

oligos_mono <- order_cells(cds = oligos_mono)

#Plots #####
plot_cells(oligos_mono,
           color_cells_by = "pseudotime",
           label_cell_groups=F,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5, show_trajectory_graph = T)

#extract pseudotime coordinates
traj.coord <- oligos_mono@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
oligos@meta.data$pseudotime <- traj.coord

#__distribution plots ####
meta_data <- oligos@meta.data %>%
  mutate(pseudo_perc = pseudotime/max(pseudotime)) %>% 
  group_by(orig.ident) %>% 
  arrange(pseudotime) %>% 
  mutate(pseudo_cum = cumsum(pseudo_perc),
         pseudo_cum_perc = pseudo_cum/sum(pseudo_cum)) %>% 
  select(orig.ident, Sex, Condition, Experiment, Batch, named.clusters, seurat_clusters, pseudotime, pseudo_cum, pseudo_cum_perc, pseudo_perc) %>% 
  group_by("orig.ident") %>% 
  mutate(names = case_when(orig.ident == "ars.female.control" ~ "Baseline\nFemale Control",
                           orig.ident == "ars.female.stress" ~ "Baseline\nFemale ARS",
                           orig.ident == "ars.male.control" ~ "Baseline\nMale Control",
                           orig.ident == "ars.male.stress" ~ "Baseline\nMale ARS",
                           orig.ident == "cms.female.control" ~ "CMS\nFemale Control",
                           orig.ident == "cms.female.stress" ~ "CMS\nFemale ARS",
                           orig.ident == "cms.male.control" ~ "CMS\nMale Control",
                           orig.ident == "cms.male.stress" ~ "CMS\nMale ARS"),
         names = factor(names, levels = c("Baseline\nFemale Control", "Baseline\nFemale ARS","Baseline\nMale Control", "Baseline\nMale ARS","CMS\nFemale Control","CMS\nFemale ARS","CMS\nMale Control", "CMS\nMale ARS")),
         Experiment_g = case_when(Experiment == "ars" ~ "Baseline",
                                  Experiment == "cms" ~ "CMS"))

ggplot(meta_data, aes(x = pseudotime, fill = names, alpha = 0.5)) +
  geom_density(stat = "density") +
  scale_fill_manual(values = col_full, aesthetics = "fill", name = "") +
  scale_alpha(guide = F) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_grid(Experiment~tools::toTitleCase(Sex)) +
  theme_classic() +
  ylab("Density") +
  xlab("Pseudotime") +
  theme(text = element_text(family = "Helvetica"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 14))

# cumulative distribution
plots <- meta_data %>% 
  left_join(col_pal_full, by = c("orig.ident" = "col_cat_full")) %>% 
  group_by(Batch) %>% 
  nest() %>% 
  mutate(plots = purrr::map(data, function(d){
    ggplot(d, aes(x = pseudotime, color = names)) +
      stat_ecdf(geom = "step", size = 1.5) +
      scale_color_manual(values = c(as.character(unique(d[d$Condition == "control",]$col_pan_full)), as.character(unique(d[d$Condition == "stress",]$col_pan_full))), aesthetics = "color", name = "") +
      scale_alpha(guide = F) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_x_continuous(expand = c(0, 0)) +
      xlab("Pseudotime") +
      ylab("Cumulative distribution") +
      # geom_vline(xintercept = 0.95) +
      theme_classic() +
      theme(panel.border = element_rect(size = 1 , color = "black", fill = NA),
            axis.line = element_blank(),
            legend.position = "none",
            text = element_text(family = "Helvetica"),
            axis.title = element_text(size = 24),
            axis.text = element_text(size = 20),
            strip.text = element_text(size = 20),
            legend.text = element_text(size = 14), 
            strip.background = element_rect(color = NA)
  )
  }))


#Stats #####
###__cumulative ####
#the function ecdf returns the function of the ecdf of the input, not the transformed data points! to have the transformed datapoints, call them within the ecdf function

#kolgomorov-smirnoff

data_stats <- meta_data %>% 
  ungroup() %>% 
  select(orig.ident, Sex, Condition, Experiment, Batch, pseudo_cum) %>%
  group_by(orig.ident, Sex, Condition, Experiment, Batch) %>% 
  nest() 

d <- data_stats$data[[1]]

meta_data %>% 
  group_by(Batch) %>% 
  nest()  -> tt

tt %>% 
  mutate(
    ks = purrr::map(data, function(c){
      broom::tidy(ks.test(x = c[c$Condition == "control",]$pseudo_perc,
                          y = c[c$Condition == "stress",]$pseudo_perc,
                          exact = NULL,
                          alternative = "two.sided"))
    })
  ) %>% 
  select(-data) %>% 
  unnest(cols = ks) %>% 
  ungroup() %>% 
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  write.table("KS_stats.txt", row.names = F)
