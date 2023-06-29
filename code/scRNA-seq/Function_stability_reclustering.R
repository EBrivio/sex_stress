# Brivio et al.
# scRNA-Seq on mouse PVN with 10x Genomics
# Samples:
# Controls vs Acute restraint (ARS) in M, F and OVX : ars.male.control or . stress
# Controls vs chronic mild stress (CMS) + ARS in M, F : cms.males.control or .stress


# Inputs:
#   dataset to cluster
#   resolution and PCs to test


subsampling_clusters_stability<- function(res, 
                                          pcs,
                                          seurat_object,
                                          iterations,
                                          cell_idents = "cell.cat",
                                          cell_cluster = "Neurons",
                                          cells_subsampled, # 0-1
                                          nfeatures = 4000, #variables features for clustering
                                          output_dir
){
  # Performs subsampling for a seurat object to calculate cluster stability 
  library(Seurat)
  library(tidyverse)
  # Clustering w subsampling! ######
  test_table <- crossing(res, pcs)
  
  
  results <- data_frame(cluster = as.integer(),
                        max = as.double(),
                        repetition = as.integer(),
                        pcs = as.double(),
                        res = as.double(),
                        size = as.double(),
                        full_size = as.double())
  
  #__iterating through res and pcs

  for (i in 1:nrow(test_table)){
    neurons1 <- seurat_object
    r = test_table$res[i]
    pcs = test_table$pcs[i]
    message("Starting test resolution: ", r, " and PCs: ", pcs)
    
    ## clustering full dataset #####
    neurons1 <- FindNeighbors(object = neurons1, dims = 1:pcs, verbose = F)
    neurons1 <- FindClusters(object = neurons1,
                             resolution = r, 
                             verbose = F)
    neurons1 <- RunUMAP(neurons1, dims = 1:pcs, verbose = F)
    # DimPlot(neurons1, label = T)
    n = ncol(seurat_object@assays$RNA)
    cells = round(n * cells_subsampled)
    
    # iterations through subsampling #####

    for (s in 1:iterations){
      
      ##selection of cells_subsampled of random cells####
      Idents(seurat_object) <- cell_idents
      selected_cells <- WhichCells(seurat_object, idents = cell_cluster, downsample = cells, seed = s)
      
      # DimPlot(neurons)
      # DimPlot(neurons, cells = selected_cells)
      test <- subset(seurat_object, cells = selected_cells)
      
      #PCA
      test <- ScaleData(object = test, vars.to.regress = NULL, verbose = F)
      test <- FindVariableFeatures(object = test, nfeatures = nfeatures, selection.method = 'vst', verbose = F)
      test <- RunPCA(object = test, pc.genes = test@var.genes, verbose = F)
      
      ## clustering subsampled dataset (test)
      test <- FindNeighbors(object = test, 
                            dims = 1:pcs, verbose = F)
      test <- FindClusters(object = test,
                           resolution = r, 
                           verbose = F)
      test <- RunUMAP(test, dims = 1:pcs, verbose = F)
      # DimPlot(test, label = T)
      
      #definition of maximum fraction of cells from each of the original clusters found in any of the clusters generated in the bootstrap
      # original dataset = neurons1
      # subsampled dataset = test
      
      #iteration across original clusters
      cl <- 0:(length(unique(neurons1$seurat_clusters))-1)
      # c = 0
      for (c in cl){
        #cluster to test
        # message(paste0("Resolution ", r,", PCs ", pcs, " - repetition ", s," out of 20 - calculating cluster ", c, " of ", max(cl)))
        original_clusters <- neurons1$seurat_clusters
        original_clusters <- original_clusters[original_clusters==c]
        full_size <- length(original_clusters) #size of the original cluster
        
        original_clusters <- original_clusters[names(original_clusters) %in% colnames(test@assays$RNA@data)]
        # size = sum(names(original_clusters) %in% colnames(test@assays$RNA@data))
        size = length(original_clusters)  #the number of cells from the original clusters that were picked in this random sampling, needed to calculate stability
        #name of all cells in cluster c
        cells_test <- names(original_clusters)
        
        #split test cluster in list of cluster with cells for each
        sub_clusters <- test$seurat_clusters
        list_sub <- sub_clusters %>% 
          as.data.frame() %>% rownames_to_column("cell") %>%
          rename("cluster" = ".") %>% 
          group_by(cluster) %>% 
          nest() %>% 
          arrange(cluster) %>% 
          mutate(count = map_dbl(data, function(d){
            sum(cells_test %in% d$cell)
          }))
        #calculate the max value across the list (I don't care which subcluster it corresponds to)
        results1 <- data_frame(cluster = c,
                               max = max(list_sub$count)/size,
                               repetition = s,
                               pcs = pcs,
                               res = r,
                               size = size, 
                               full_size = full_size)
        
        results <- results %>% 
          bind_rows(results1)
        
      } #end loop clusters
      message("Resolution ",r,", PCs ", pcs, " - repetition ", s," - FINISHED ALL CLUSTERS")
      
    } #end loop iterations
    message("Resolution ", r,", PCs ", pcs, " - FINISHED ALL REPETITIONS for ", i, " out of ", nrow(test_table))
    message(Sys.time())
  } #end loop resolution
  
  
  
  saveRDS(object = results, file = paste0(output_dir, "subsampling_results_res_pcs.Rds"))
  return(results)
  
}


