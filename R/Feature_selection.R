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
                                 group = 'sample', group1 = NULL, group2 = NULL,
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
    df_neighbor$batch = seuInt[group]
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
    rownames(sampletable) = sampletable[group]
    label = sampletable[colnames(res), ][group]
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
interaction_gen_continuous = function(neighbor_mat,
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
Feature_selection = function(neighbor_mat, neighbor_mat_list, prop_df, markers, type,
                             n_cluster = 8, zero_rm = T, ...){
  if(type == 'grouped'){
    interaction_gen = interaction_gen_group
  }else if(type == 'continuous'){
    interaction_gen = interaction_gen_continuous
  }else{
    warning("Type must be either 'grouped' or 'continuous'")
  }
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
