### Figure 3
seurat_object <- readRDS('combine_dataset.rds')

gene_count_matrix <- as.matrix(seurat_object@assays$RNA@counts)
### Fig. 3 dot plot
select_genes <- c("mro",'galT','galK','glgX','glgX','treS','susB','pgcA','glkA','pckA','pyk',
                  'mdh','ttdA','fumB','ldhA','por','cutC','cat1','dhaT',
                  'pgi','pfkA','fbp','fba','fbaA','fda','tpiA','gap','gapA','pgk','gpmI','gpmA','eno',
                  'zwf','gnd',"rpiB",'tkt','tal','prs','gltA','acnD','icd','frdA')

bubble.df <- as.data.frame(t(gene_count_matrix))

bubble.df <- bubble.df[,select_genes]

bubble.df <- as.data.frame(bubble.df)

rownames(bubble.df)==rownames(seurat_object@meta.data) -> id

table(id)

bubble.df_avg_byspecies <- aggregate(bubble.df[,1:ncol(bubble.df)],list(seurat_object$species_info),mean)

rownames(bubble.df_avg_byspecies) <- bubble.df_avg_byspecies$Group.1

bubble.df_avg_byspecies <- bubble.df_avg_byspecies[,-1]

bubble.df_avg_byspecies <- as.data.frame(t(bubble.df_avg_byspecies))


pheatmap(bubble.df_avg_byspecies, 
         
         color = pal,
         
         border_color = NA,#scale = 'row',
         
         #labels_col="",
         
         cluster_cols = F,show_colnames = T,cluster_rows = F,treeheight_row = 0,
         
         #gaps_col = c(7,14),
         #cutree_rows = 2,
         fontsize=10, fontsize_row=10,fontsize_col=10) #自定义颜色
