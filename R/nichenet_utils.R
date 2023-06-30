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
#'
Generate_nichenet = function(seuInt, sender, receiver, label = 'label', group1, group2){
  new_obj = CreateSeuratObject(counts = seuInt[['count']])
  new_obj@meta.data = seuInt@meta.data
  Idents(new_obj) = 'cluster'
  nichenet_output = nichenet_seuratobj_aggregate(
    seurat_obj = new_obj,
    receiver = 1,
    condition_colname = label, condition_oi = "PN", condition_reference = "HC",
    sender = 2:7,
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network,
    weighted_networks = weighted_networks, organism = "human",
    expression_pct = 0, lfc_cutoff = 0, geneset = "DE", filter_top_ligands = TRUE ,
    top_n_ligands = 100, top_n_targets = 500, cutoff_visualization = 0)

  nichenet_output$top_ligands
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
#'
plot_nichenet = function(markers, p_list, plot = T){
  p_list = p_list %>% filter(cluster == 'all') %>% filter(neighbor_cluster == 'all')
  senders = intersect(rownames(ligand_target_matrix), markers)
  receivers = intersect(colnames(ligand_target_matrix), markers)

  mat_niche = ligand_target_matrix[senders, receivers]
  mat_p = mat_niche
  for(i in 1:length(senders)){
    sender = senders[i]
    p_list_sub = p_list[p_list[,2]==sender, ]
    rownames(p_list_sub) = p_list_sub[,1]
    p_list_sub = p_list_sub[receivers, ]
    mat_p[sender, ] = p_list_sub$p
  }

  df_niche = data.frame(melt(mat_niche))
  df_p = data.frame(melt(mat_p))

  df = as.data.frame(cbind(log(df_niche$value), log(df_p$value)))
  colnames(df) = c('Nichenet', 'p')

  ggplot(df, aes(x = Nichenet, y = p)) +
    geom_point() + geom_smooth(method='lm', formula= y ~ x, se = F) +
    ggtitle(paste0('Correlation', cor(df$Nichenet, df$p)))
  ggsave(paste0('.//figures//', 'correlation', '.png'),
         width = 20, height = 18, units = "cm")

}
