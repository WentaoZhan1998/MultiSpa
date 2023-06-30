# MultiSpa

```{r}
source('/dcs05/hongkai/data/wzhan/Spatial_visium/pipeline.R')
read_path = "/dcs05/hongkai/data/wzhan/Spatial_visium/Liver/data/"
save_path = "/dcs05/hongkai/data/wzhan/Spatial_visium/Liver/"
files = list.files(read_path)
impute = T

obj_list = create_obj_list(read_path, files, impute = impute)
## Create an object that has
## assays 'count', 'PRE_CAST';
## reductions 'PRECAST', 'position', 'umap', 'pca', 'harmony'(optional)
## clusters 'cluster', 'seurat_clusters'(optional)
seuInt = PRECAST_pipeline(obj_list, premin.spots = 0, premin.features = 0,
                          postmin.spots = 0, postmin.features = 0)
seuInt = Create_sample(seuInt, files)
saveRDS(seuInt, paste0(save_path, 'seuInt', ifelse(impute, '_impute', ''), '.rds'))

seuInt = readRDS(paste0(save_path,'seuInt.rds'))
seuInt_imput = readRDS(paste0(save_path,'seuInt_impute.rds'))

variable_genes = intersect(rownames(seuInt_imput), rownames(seuInt))
pairs = lr_load("fantom5",NULL,'human',variable_genes)
markers = unlist(lr_load("fantom5",NULL,'human',variable_genes))
names(markers) = NULL
markers = unique(markers)
saveRDS(list(markers = markers, pairs = pairs), paste0(save_path, 'markers_shared.rds'))

Gene_mat_generate(seuInt, markers, k=20,
                  save.path = save_path, name = "neighbor_gene_NICHES_k20_shared")
Gene_mat_generate(seuInt_imput, markers, k=20,
                  save.path = save_path, name = "neighbor_gene_imput_NICHES_k20_shared")

impute = F
seuInt = readRDS(paste0(save_path,'seuInt', ifelse(impute, '_impute', ''), '.rds'))
neighbor_mat = readRDS(paste0(save_path,'neighbor_gene',
                              ifelse(impute, '_imput', ''), '_NICHES_k20_shared.rds'))
neighbor_mat_list = readRDS(paste0(save_path,'neighbor_gene',
                                   ifelse(impute, '_imput', ''), '_NICHES_k20_shared_list.rds'))
markers = readRDS(paste0(save_path, 'markers_shared.rds'))$markers
pairs = readRDS(paste0(save_path, 'markers_shared.rds'))$pairs
expr_mat = t(GetAssayData(seuInt, slot = "data", assay = 'count')[markers, ])
colnames(expr_mat) = toupper(colnames(expr_mat))

##### Neighbor construction #####
prop_df = neighbor_mat_gen(seuInt, k = 20, save.path = save_path,
                           name = paste0('neighbor_cluster_k20', ifelse(impute, '_impute', '')))$prop_df
prop_df = readRDS(paste0(save_path, 'neighbor_cluster_k20',
                         ifelse(impute, '_impute', ''), '.rds'))
n_cluster = length(unique(prop_df$cluster))

sampletable = data.frame(batch = 1:length(files), sample = files)
sampletable$label = rep(c('L', "N", "T"), 4)
sampletable$label_num = rep(c(1, 0, 2), 4)
sampletable$batch = rep(1:4, each = 3)
seuInt$label = left_join(data.frame(sample = seuInt$sample), sampletable, by = 'sample')$label
seuInt$label_num = left_join(data.frame(sample = seuInt$sample), sampletable, by = 'sample')$label_num
seuInt$batch = left_join(data.frame(sample = seuInt$sample), sampletable, by = 'sample')$batch

##### Top features selection #####
zero_rm = F
p_top = Feature_selection(neighbor_mat, neighbor_mat_list, prop_df, markers = colnames(neighbor_mat),
                          interaction_gen = interaction_gen_group, zero_rm = zero_rm ,
                          seuInt = seuInt, expr_mat = expr_mat,
                          sampletable = sampletable,
                          group1 = c('HCC-1N', 'HCC-2N', 'HCC-3N', 'HCC-4N'),
                          group2 = c('HCC-1T', 'HCC-2T', 'HCC-3T', 'HCC-4T'))
p_top = arrange(p_top, p.adj)
head(p_top)
saveRDS(p_top, paste0(save_path, 'p_top_raw',
                      ifelse(impute, '_impute', ''),
                      ifelse(zero_rm, '_zerorm', '_zerokept'),
                      '_k20_shared.rds'))

p_top = Feature_selection(neighbor_mat, neighbor_mat_list, prop_df, markers = colnames(neighbor_mat),
                          interaction_gen = interaction_gen_time, zero_rm = zero_rm ,
                          seuInt = seuInt, expr_mat = expr_mat,
                          sampletable = sampletable,
                          group = 'label_num')
p_top = arrange(p_top, p.adj)
head(p_top)
saveRDS(p_top, paste0(save_path, 'p_time_raw',
                      ifelse(impute, '_impute', ''), '_zerokept_k20_shared.rds'))

##### Top features plots #####
p_top = readRDS(paste0(save_path, 'p_time_raw',
                       ifelse(impute, '_impute', ''),
                       ifelse(zero_rm, '_zerorm', '_zerokept'),
                       '_k20_shared.rds'))
for(i in 1:length(pairs[[1]])){
  print(i)
  id1 = p_top[,1] == pairs[[1]][i]
  id2 = p_top[,2] == pairs[[2]][i]
  #temp = p_top[id1&id2,] %>% filter(cluster=='all') %>% filter(neighbor_cluster == 'all')
  temp = p_top[id1&id2,]
  temp = temp %>% filter(cluster == 'all') %>% filter(neighbor_cluster == 'all')
  if(nrow(temp) > 0){
    for(j in 1:nrow(temp)){
      plot_gene_by_gene(pairs[[1]][i], pairs[[2]][i], temp[j,5], temp[j,6], p_adj = temp$p[j],
                        expr_mat, neighbor_mat, neighbor_mat_list,
                        cell_by = seuInt$label, smooth_by = seuInt$label, facet_by = seuInt$batch, zero_rm = zero_rm)
    }
  }
}

for(sample in unique(seuInt$sample)){
  #for(sample in c('SL1-A1', 'SL4-D1')){
  test = Load10X_Spatial(
    paste0(read_path, sample),
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial"
  )
  test$sample = rep(sample, ncol(test))
  test = subset(test, subset = nCount_Spatial > quantile(test@meta.data$nCount_Spatial, 0.05))
  test = subset(test, subset = nFeature_Spatial > quantile(test@meta.data$nFeature_Spatial, 0.05))
  for(gene in c('VWF', 'SIRPA')){
    print(gene)
    SpatialFeaturePlot(test, features = gene, pt.size.factor = 2)
    ggsave(paste0(save_path, 'figures/', sample, '_', gene, '.png'))
  }
}

process_fun = function(x){
  lm(neighbor_gene ~ gene, data = x)$coef[2]
}

res = sapply(split(df_plt, df_plt$group_smooth), process_fun)
```
