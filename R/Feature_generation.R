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
pseudobulk_gen = function(seuInt, markers){
  obj_list = SplitObject(seuInt, split.by = "sample")
  pseudobulk = t(
    do.call(rbind,
            lapply(obj_list, function(obj){
              melt(AggregateExpression(obj, features = markers, assay = 'count')[[1]])$value
            }))
  )
  colnames(pseudobulk) = names(obj_list)
  rownames(pseudobulk) = paste(
    melt(AggregateExpression(seuInt, features = markers, assay = 'count')[[1]])[,1],
    melt(AggregateExpression(seuInt, features = markers, assay = 'count')[[1]])[,2],
    sep = '-')
  return(pseudobulk)
}




