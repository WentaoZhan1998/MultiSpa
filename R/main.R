##### Use any python environment proper for the loaded R package

my_tolower = function(x){
  x_split = tolower(unlist(strsplit(x, split = '')))
  x_split[1] = toupper(x_split[1])
  return(paste(x_split, collapse = ''))
}

#' Sequence optimize
#'
#' Do sequence optimize
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#' @param model the model to be used. The default is set as 'linear' to use the linear model.
#' @param option a list specifies the options for model solving, currently supporting 'mle_solver' only.
#'
#' @return A list including the information of the fitted model
#'
#' @export
create_obj_list = function(path, files = NULL, impute = T, harmony = T){
  obj_list = list()
  names = vector()
  if(is.null(files)){
    files = list.files(path)
  }
  l = length(files)

  for(i in 1:l){
    file = files[i]
    message(file)
    test = Load10X_Spatial(
      paste0(path, file),
      filename = "filtered_feature_bc_matrix.h5",
      assay = "Spatial"
    )
    test$sample = rep(file, ncol(test))
    test = subset(test, subset = nCount_Spatial > quantile(test@meta.data$nCount_Spatial, 0.05))
    test = subset(test, subset = nFeature_Spatial > quantile(test@meta.data$nFeature_Spatial, 0.05))
    test$col = test@images$slice1@coordinates$imagerow
    test$row = test@images$slice1@coordinates$imagecol
    if(impute){
      test_denoise <- magic(t(as.matrix(GetAssayData(test))), genes="all_genes", t=6)
      #test_denoise <- saver(t(readRDS('test.rds')))
      #test_denoise <- VIPER(as.matrix(GetAssayData(test)), num = 10, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1,
      #             report = FALSE, outdir = NULL, prefix = NULL)
      test[['RNA']] = CreateAssayObject(t(test_denoise$result))
      test[['Spatial']] = CreateAssayObject(t(test_denoise$result))
    }
    else{
      test[['RNA']] = CreateAssayObject(GetAssayData(test))
      test[['Spatial']] = CreateAssayObject(GetAssayData(test))
    }
    obj_list = append(obj_list, test)
  }
  return(obj_list)
}

#### Precesss with PRECAST
#' Sequence optimize
#'
#' Do sequence optimize
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#' @param model the model to be used. The default is set as 'linear' to use the linear model.
#' @param option a list specifies the options for model solving, currently supporting 'mle_solver' only.
#'
#' @return A list including the information of the fitted model
#'
#' @export
PRECAST_pipeline = function(obj_list, k = 8, harmony = T, ...){
  set.seed(2022)
  PRECASTObj <- CreatePRECASTObject(obj_list, ...)
  seuCount = Merge_Seurat_List(PRECASTObj@seulist, add.cell.ids = NULL,
    merge.data = TRUE, project = "SeuratProject")
  ### Add the model setting
  ## check the number of genes/features after filtering step
  ## Add adjacency matrix list for a PRECASTObj object to prepare for PRECAST model fitting.
  PRECASTObj <-  AddAdjList(PRECASTObj, platform = "Visium")
  ## Add a model setting in advance for a PRECASTObj object. verbose =TRUE helps outputing the information in the algorithm.
  PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal=FALSE, verbose=TRUE, int.model=NULL)
  PRECASTObj <- PRECAST(PRECASTObj, K=k)
  PRECASTObj <- selectModel(PRECASTObj)
  seuInt <- IntegrateSpaData(PRECASTObj, species='Human')
  seuInt <- ScaleData(seuInt, verbose = FALSE)
  seuCount = RenameCells(seuCount, new.names = colnames(seuInt))
  count <- GetAssayData(object =  seuCount[['Spatial']])
  seuInt[["count"]] <- CreateAssayObject(data = count)
  DefaultAssay(seuInt) = 'count'
  seuInt = FindVariableFeatures(seuInt, nfeatures = 2000)
  seuInt <- ScaleData(seuInt)
  seuInt <- RunPCA(seuInt, features = VariableFeatures(object = seuInt))
  seuInt <- FindNeighbors(seuInt, reduction = 'PRECAST', graph.name = c('PRECAST_nn', 'PRECAST_snn'))
  if(harmony){
    seuInt <- RunHarmony(seuInt, 'batch', reduction = 'pca', max.iter.harmony = 100)
    seuInt <- FindNeighbors(seuInt, reduction = 'harmony', graph.name = c('harmony_nn', 'harmony_snn'))
    seuInt <- FindClusters(seuInt, graph.name = 'harmony_snn', resolution = 0.4) #Can work independently
  }
  Idents(seuInt) = 'cluster'

  return(seuInt)
}

#### Create sample label for seurat object
#' Sequence optimize
#'
#' Do sequence optimize
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#' @param model the model to be used. The default is set as 'linear' to use the linear model.
#' @param option a list specifies the options for model solving, currently supporting 'mle_solver' only.
#'
#' @return A list including the information of the fitted model
#'
#' @export
Create_sample = function(seuInt, files){
  sampletable = data.frame(batch = 1:length(files), sample = files)
  seuInt$sample = left_join(data.frame(batch = as.integer(seuInt$batch)), sampletable, by = 'batch')$sample
  return(seuInt)
}

#### Create spatial plots, with cells labeled by the cluster IDs.
#' Sequence optimize
#'
#' Do sequence optimize
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#' @param model the model to be used. The default is set as 'linear' to use the linear model.
#' @param option a list specifies the options for model solving, currently supporting 'mle_solver' only.
#'
#' @return A list including the information of the fitted model
#'
#' @export
Colorplot = function(seuInt, save.path, name = "Int_cluster_", cluster = NULL,
                     cols_cluster = c( "#FD7446" ,"#709AE1", "#31A354","#9EDAE5",
                                       "#DE9ED6" ,"#BCBD22", "#CE6DBD" ,"#DADAEB" ,
                                       "yellow", "#FF9896","#91D1C2", "#C7E9C0" ,"#6B6ECF", "#7B4173" )){
  if(is.null(cluster)){Idents(seuInt) = 'cluster'}
  pList = SpaPlot(seuInt, batch=NULL,point_size=1.5, cols=cols_cluster, combine=F)
  files = unique(seuInt$sample)
  l = length(files)
  j = 1
  for(i in c(1:l)){
    file = files[i]
    p = pList[[j]] + ggtitle(label = file)
    print(p)
    ggsave(paste0(save.path, name, file, ".png"), width = 12, height = 10, units = "cm")
    j = j+1
  }
}

#### Plot with PRECAST's cluster label
Reduce = function(seuInt, reduction = 'PRECAST',
                  cluster = 'cluster', sample = 'sample',
                  plot = T, save.path, name = "_Umap"){
  ##### PRECAST's cluster label #####
  seuInt <- ScaleData(seuInt, verbose = FALSE)
  seuInt <- FindNeighbors(seuInt, reduction = reduction, dims = 1:15)
  seuInt <- RunUMAP(seuInt, reduction = reduction, dims = 1:15, min.dist = 0.3)
  if(plot){
    DimPlot(seuInt, reduction = "umap", group.by = cluster, label = TRUE, repel = TRUE)
    ggsave(paste0(save.path, name,".png"), width = 14, height = 10, units = "cm")

    DimPlot(seuInt, reduction = "umap", group.by = sample, label = TRUE, repel = TRUE)
    ggsave(paste0(save.path, reduction, name,"_sample.png"), width = 14, height = 10, units = "cm")
  }
}

#### Generate and save the neighboring expression matrix
Gene_mat_generate_simple = function(seuInt, markers, k=20){
  neighbor_mat = matrix(nrow = ncol(seuInt), ncol = length(markers))
  colnames(neighbor_mat) = markers
  for(sample0 in unique(seuInt$batch)){
    print(sample0)
    id = seuInt$batch == sample0
    DefaultAssay(seuInt) <- "count"
    seuInt_temp = subset(seuInt, batch == sample0)
    temp = GetAssayData(seuInt_temp)
    assay_temp = temp[toupper(rownames(temp)) %in% markers,]
    seuInt_temp <- FindNeighbors(seuInt_temp, reduction = 'position', dims = 1:2,  k.param=k,
                                 graph.name = c('spatial_nn', 'spatial_snn'))
    neighbor_mat_temp = apply(seuInt_temp@graphs$spatial_nn, 1, function(x){
      rowMeans(assay_temp[,x==1])
    })
    neighbor_mat[id, ] = t(neighbor_mat_temp)
  }

  return(neighbor_mat)
}

#' Sequence optimize
#'
#' Do sequence optimize
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#' @param model the model to be used. The default is set as 'linear' to use the linear model.
#' @param option a list specifies the options for model solving, currently supporting 'mle_solver' only.
#'
#' @return A list including the information of the fitted model
#'
#' @export
#'
Gene_mat_generate = function(seuInt, markers, k=20, save.path, name = NULL){
  if(is.null(name)){
    name = paste0("neighbor_gene", k)
  }
  neighbor_mat = matrix(nrow = ncol(seuInt), ncol = length(markers))
  colnames(neighbor_mat) = markers
  neighbor_mat_list = list()
  for(i in unique(seuInt$cluster)){
    neighbor_mat_list[[i]] = matrix(nrow = ncol(seuInt), ncol = length(markers))
    colnames(neighbor_mat_list[[i]]) = markers
  }
  for(sample0 in unique(seuInt$batch)){
    print(sample0)
    id = seuInt$batch == sample0
    DefaultAssay(seuInt) <- "count"
    seuInt_temp = subset(seuInt, batch == sample0)
    temp = GetAssayData(seuInt_temp)
    assay_temp = temp[toupper(rownames(temp)) %in% markers,]
    seuInt_temp <- FindNeighbors(seuInt_temp, reduction = 'position', dims = 1:2,  k.param=k,
                                 graph.name = c('spatial_nn', 'spatial_snn'))
    neighbor_mat_temp = apply(seuInt_temp@graphs$spatial_nn, 1, function(x){
      rowMeans(assay_temp[,x==1])
    })
    neighbor_mat[id, ] = t(neighbor_mat_temp)
    for(i in unique(seuInt$cluster)){
      print(paste0('cluster', as.character(i)))
      neighbor_mat_temp = apply(seuInt_temp@graphs$spatial_nn, 1, function(x){
        id_temp = (x==1 & seuInt_temp$cluster==i)
        if(sum(id_temp) == 0){
          return(rep(0, length(markers)))
        }
        else if(sum(id_temp) == 1){
          return(as.vector(assay_temp[,id_temp]))
        }
        else{
          return(rowMeans(assay_temp[,id_temp]))
        }
      })
      neighbor_mat_list[[i]][id, ] = t(neighbor_mat_temp)
    }
  }

  saveRDS(neighbor_mat, paste0(save.path, name, '.rds'))
  saveRDS(neighbor_mat_list, paste0(save.path, name, '_list.rds'))
}

#### Generate and save the neighboring proporation matrix
#' Sequence optimize
#'
#' Do sequence optimize
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#' @param model the model to be used. The default is set as 'linear' to use the linear model.
#' @param option a list specifies the options for model solving, currently supporting 'mle_solver' only.
#'
#' @return A list including the information of the fitted model
#'
#' @export
neighbor_mat_gen = function(seuInt, k = 20, save.path, name = NULL){
  if(is.null(name)){
    name = paste0("neighbor_cluster", k)
  }
  seuInt$cluster = as.integer(seuInt$cluster)
  neighbor_mat = matrix(nrow = ncol(seuInt), ncol = k)
  for(sample0 in unique(seuInt$batch)){
    print(sample0)
    id = seuInt$batch == sample0
    seuInt_temp = subset(seuInt, batch == sample0)
    seuInt_temp <- FindNeighbors(seuInt_temp, reduction = 'position', dims = 1:2,  k.param=k,
                                 graph.name = c('spatial_nn', 'spatial_snn'))
    neighbor_mat_temp = apply(seuInt_temp@graphs$spatial_nn, 1, function(x){
      seuInt_temp$cluster[x==1]
    })
    neighbor_mat[id, ] = t(neighbor_mat_temp)
  }
  prop_mat = matrix(nrow = ncol(seuInt), ncol = length(unique(seuInt$cluster)))
  for(i in 1:length(unique(seuInt$cluster))){
    cluster = i
    print(cluster)
    prop_mat[,i] = rowSums(neighbor_mat == cluster)
  }
  prop_df = data.frame(prop_mat)
  rm(prop_mat)
  colnames(prop_df) = 1:length(unique(seuInt$cluster))
  prop_df$cluster = seuInt$cluster
  prop_df$batch = seuInt$batch

  res = list()
  res$neighbor_mat = neighbor_mat
  res$prop_df = prop_df
  saveRDS(res, paste0(save.path, name, '.rds'))
  return(res)
}

#' Sequence optimize
#'
#' Do sequence optimize
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#' @param model the model to be used. The default is set as 'linear' to use the linear model.
#' @param option a list specifies the options for model solving, currently supporting 'mle_solver' only.
#'
#' @return A list including the information of the fitted model
#'
#' @export
#'
interaction_gen_group = function(neighbor_mat,
                           seuInt, expr_mat, sampletable,
                           process_fun,
                           id = NULL, id_temp_fun = NULL,
                           markers = NULL, markers2 = NULL,
                           group1 = NULL, group2 = NULL,
                           outnames = c('variable', 'neighbor_variables'),
                           expr_needed = F, zero_rm = F){
  if(is.null(markers)){
    markers = colnames(neighbor_mat)
  }
  if(is.null(markers2)){
    markers2 = markers
  }
  if(is.null(id)){
    id = rep(T, nrow(neighbor_mat))
  }
  if(is.null(id_temp_fun)){
    id_temp_fun = function(i){rep(T, nrow(neighbor_mat))}
  }
  p_mat_correlation = matrix(nrow  = length(markers), ncol = length(markers2))
  rownames(p_mat_correlation) = markers
  colnames(p_mat_correlation) = markers2
  for(i in 1:length(markers)){
    id_temp = id_temp_fun(i)
    df_neighbor = data.frame(neighbor_mat[,markers2])
    colnames(df_neighbor) = markers2
    df_neighbor$cluster = seuInt$cluster
    df_neighbor$batch = seuInt$sample
    if(expr_needed == T){
      df_neighbor$expr = expr_mat[,i]
      if(zero_rm == T){
        df_neighbor = filter(df_neighbor, expr>0)
      }
    }
    df_neighbor = df_neighbor[id & id_temp, ]
    if(sum(complete.cases(df_neighbor))<10){next}
    df_neighbor = df_neighbor[complete.cases(df_neighbor),]
    res = sapply(split(df_neighbor, df_neighbor$batch), process_fun)
    colnames(res) = names(split(df_neighbor, df_neighbor$batch))
    res[is.na(res)] = 0
    rownames(sampletable) = sampletable$sample
    label = sampletable[colnames(res), ]$sample
    if(sum(label %in% group1)<=1 |
       sum(label %in% group2)<=1){next}
    p_mat_correlation[i, ] = apply(res, 1, function(x){
      res = t.test(x[label %in% group1],
                   x[label %in% group2],
                   alternative = "two.sided", var.equal = FALSE)
      res$p.value
    })
  }
  p_df_correlation = melt(p_mat_correlation)
  colnames(p_df_correlation) = c(outnames, 'p')
  p_df_correlation$p.adj = p.adjust(p_df_correlation$p, method = "hochberg")
  p_df_correlation = arrange(p_df_correlation, p.adj)
  return(p_df_correlation)
}

overall_p <- function(my_model) {
  f <- summary(my_model)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

#' Sequence optimize
#'
#' Do sequence optimize
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#' @param model the model to be used. The default is set as 'linear' to use the linear model.
#' @param option a list specifies the options for model solving, currently supporting 'mle_solver' only.
#'
#' @return A list including the information of the fitted model
#'
#' @export
interaction_gen_time = function(neighbor_mat,
                                 seuInt, expr_mat, sampletable,
                                 process_fun,
                                 id = NULL, id_temp_fun = NULL,
                                 markers = NULL, markers2 = NULL,
                                 group = NULL,
                                 outnames = c('variable', 'neighbor_variables'),
                                 expr_needed = F, zero_rm = F){
  if(is.null(markers)){
    markers = colnames(neighbor_mat)
  }
  if(is.null(markers2)){
    markers2 = markers
  }
  if(is.null(id)){
    id = rep(T, nrow(neighbor_mat))
  }
  if(is.null(id_temp_fun)){
    id_temp_fun = function(i){rep(T, nrow(neighbor_mat))}
  }
  p_mat_correlation = matrix(nrow  = length(markers), ncol = length(markers2))
  rownames(p_mat_correlation) = markers
  colnames(p_mat_correlation) = markers2
  for(i in 1:length(markers)){
    id_temp = id_temp_fun(i)
    df_neighbor = data.frame(neighbor_mat[,markers2])
    colnames(df_neighbor) = markers2
    df_neighbor$cluster = seuInt$cluster
    df_neighbor$batch = seuInt$sample
    if(expr_needed == T){
      df_neighbor$expr = expr_mat[,i]
      if(zero_rm == T){
        df_neighbor = filter(df_neighbor, expr>0)
      }
    }
    df_neighbor = df_neighbor[id & id_temp, ]
    if(sum(complete.cases(df_neighbor))<10){next}
    df_neighbor = df_neighbor[complete.cases(df_neighbor),]
    res = sapply(split(df_neighbor, df_neighbor$batch), process_fun)
    colnames(res) = names(split(df_neighbor, df_neighbor$batch))
    res[is.na(res)] = 0
    rownames(sampletable) = sampletable$sample
    time = sampletable[colnames(res), ][[group]]
    if(length(unique(time))==1){next}
    p_mat_correlation[i, ] = apply(res, 1, function(x){
      model = lm(x ~ time)
      overall_p(model)
    })
  }
  p_df_correlation = melt(p_mat_correlation)
  colnames(p_df_correlation) = c(outnames, 'p')
  p_df_correlation$p.adj = p.adjust(p_df_correlation$p, method = "hochberg")
  p_df_correlation = arrange(p_df_correlation, p.adj)
  return(p_df_correlation)
}

#' Sequence optimize
#'
#' Do sequence optimize
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#' @param model the model to be used. The default is set as 'linear' to use the linear model.
#' @param option a list specifies the options for model solving, currently supporting 'mle_solver' only.
#'
#' @return A list including the information of the fitted model
#'
#' @export
Feature_selection = function(neighbor_mat, neighbor_mat_list, prop_df, markers, interaction_gen,
                             n_cluster = 8, zero_rm = T, ...){
  p_top = data.frame(matrix(nrow = 0, ncol = 6))
  if(T){
    p_df_temp = interaction_gen(prop_df[, 1:(ncol(prop_df)-2)],
                                process_fun = function(x){
                                  colMeans(x[ , !names(x) %in% c('expr', 'batch', 'cluster')])
                                },
                                id = NULL, id_temp_fun = function(i){seuInt$cluster == i},
                                markers = as.character(1:n_cluster), markers2 = NULL,
                                zero_rm = zero_rm, expr_needed = F, ...)
    p_df_temp$cluster = 'all'
    p_df_temp$neighbor_cluster = 'all'
    p_top = rbind(p_top, p_df_temp)
  }
  for(i in c(names(neighbor_mat_list), 'all')){
    if(i == 'all'){
      neighbor_temp = neighbor_mat
      id_temp1 = NULL
    }
    else{
      neighbor_temp = neighbor_mat_list[[i]]
      id_temp1 = seuInt$cluster == i
    }
    if(T){
      p_df_temp = interaction_gen(neighbor_temp,
                                  process_fun = function(x){
                                    colMeans(x[ , !names(x) %in% c('expr', 'batch', 'cluster')])
                                  },
                                  id = NULL, id_temp_fun = function(i){seuInt$cluster == i},
                                  markers = as.character(1:n_cluster), markers2 = markers,
                                  zero_rm = zero_rm, expr_needed = F, ...)
      p_df_temp$cluster = 'all'
      p_df_temp$neighbor_cluster = i
      p_top = rbind(p_top, p_df_temp)
      p_df_temp = interaction_gen(prop_df[,1:(ncol(prop_df)-2)],
                                  process_fun = function(x){
                                    cor(x$expr, x[ , !names(x) %in% c('expr', 'batch', 'cluster')])
                                  },
                                  id = id_temp1, id_temp_fun = NULL,
                                  markers = markers, markers2 = as.character(1:n_cluster),
                                  zero_rm = zero_rm, expr_needed = T, ...)
      p_df_temp$cluster = i
      p_df_temp$neighbor_cluster = 'all'
      p_top = rbind(p_top, p_df_temp)
    }
    for(j in c(names(neighbor_mat_list), 'all')){
      if(j == 'all'){
        id_temp2 = NULL
      }
      else{
        id_temp2 = seuInt$cluster == j
      }
      p_df_temp = interaction_gen(neighbor_temp,
                                  process_fun = function(x){
                                    cor(x$expr, x[ , !names(x) %in% c('expr', 'batch', 'cluster')])
                                  },
                                  id = id_temp2, id_temp_fun = NULL,
                                  markers = markers, markers2 = NULL,
                                  zero_rm = zero_rm, expr_needed = T, ...)
      p_df_temp$cluster = j
      p_df_temp$neighbor_cluster = i
      p_top = rbind(p_top, p_df_temp)
      print(paste0(i, j))
    }
  }
  p_top = arrange(p_top, p.adj)
  return(p_top)
  #saveRDS(p_top, 'p_top_raw_zerokept_k20_shared.rds')
}

#' Sequence optimize
#'
#' Do sequence optimize
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#' @param model the model to be used. The default is set as 'linear' to use the linear model.
#' @param option a list specifies the options for model solving, currently supporting 'mle_solver' only.
#'
#' @return A list including the information of the fitted model
#'
#' @export
plot_cluster_gen = function(cluster, gen, p_adj, seuInt){
  neighbor_temp = neighbor_mat
  df_plt = data.frame(gene = neighbor_mat[,gen], label = seuInt$label)[seuInt$cluster == cluster, ]
  ggplot(data = df_plt, aes(fill=label, y=gene, x=label)) +
    geom_violin(position="dodge") +
    xlab("") + ylab("Tip (%)") +
    ggtitle(paste0('gene', gen, ' in celltype ', as.character(cluster), 's neighbor ,
                   , p_adj=', as.character(p_adj)))
  ggsave(paste0('.//figures//', cluster, ' vs ', gen, '.png'),
         width = 20, height = 18, units = "cm")
}

#' Sequence optimize
#'
#' Do sequence optimize
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#' @param model the model to be used. The default is set as 'linear' to use the linear model.
#' @param option a list specifies the options for model solving, currently supporting 'mle_solver' only.
#'
#' @return A list including the information of the fitted model
#'
#' @export
plot_gen_cluster = function(gen, cluster, p_adj, seuInt, expr_mat, zero_rm = T){
  neighbor_temp = prop_df
  df_plt = data.frame(gene = expr_mat[, gen], cluster = prop_df[ ,cluster], label = seuInt$label, smooth = seuInt$sample)
  if(zero_rm){
    df_plt = df_plt[df_plt$gene>0, ]
  }
  ggplot(df_plt, aes(x = cluster, y = gene, color = label)) +
    geom_point() +
    geom_smooth(aes(group = smooth), method='lm', formula= y ~ x, se = F) +
    ggtitle(paste0('gene ', gen,
                   ' vs ', 'celltype ', as.character(cluster),
                   ' proportion in neighbor , p_adj=', as.character(p_adj)))
  ggsave(paste0('.//figures//', gen, ' vs ', cluster, ifelse(zero_rm, '_zerorm', 'zero_kept'), '.png'),
         width = 20, height = 18, units = "cm")
}

#' Sequence optimize
#'
#' Do sequence optimize
#'
#' @details
#'
#' @param design matrix, the design matrix X for the predictors.
#' @param outcome vector, the output Y for the response.
#' @param model the model to be used. The default is set as 'linear' to use the linear model.
#' @param option a list specifies the options for model solving, currently supporting 'mle_solver' only.
#'
#' @return A list including the information of the fitted model
#'
#' @export
plot_gene_by_gene = function(gene1, gene2, i, j, p_adj,
                             expr_mat, neighbor_mat, neighbor_mat_list,
                             cell_by = seuInt$label, smooth_by = seuInt$sample, facet_by = seuInt$batch,
                             zero_rm = F){
  if(i == 'all'){
    neighbor_temp = neighbor_mat
  }
  else{
    neighbor_temp = neighbor_mat_list[[i]]
  }
  if(j == 'all'){
    id_temp = rep(T, nrow(seuInt))
  }
  else{
    id_temp = seuInt$cluster == j
  }
  df_plt = data.frame(gene = expr_mat[,gene1], neighbor_gene = neighbor_temp[,gene2],
                      group_cell = cell_by, group_smooth = smooth_by, group_facet = facet_by)[id_temp, ]
  if(zero_rm == T){
    df_plt = filter(df_plt, gene>0)
  }
  k = length(unique(df_plt$group_cell))
  ggplot(df_plt, aes(x = neighbor_gene, y = gene)) +
    geom_point(aes(color = as.factor(group_cell)), alpha = 0.5, size = 1) +
    geom_smooth(aes(group = group_smooth, color = as.factor(group_cell)), method='lm', formula= y ~ x, se = F) +
    #scale_color_manual(values = hue_pal()(2*k)[2*(1:k)-1]) +
    #facet_wrap(group_cell~group_facet, ncol = 2) +
    ggtitle(paste0('cluster ', as.character(i), ' ', gene1,
                   ' vs ', 'cluster ', as.character(j), ' ', gene2,
                   ' in neighbor , p=', as.character(p_adj)))
  ggsave(paste0('.//figures//', gene1, ' vs ', gene2, ' ', as.character(i), as.character(j), '.png'),
         width = 20, height = 18, units = "cm")
}


