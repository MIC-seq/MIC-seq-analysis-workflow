### Figure 6

## Fig. 6A
DimPlot(mega_seurat, reduction = "umap",label.size = 6,label = F,cols = cell_type_cols,
        
        group.by = 'seurat_clusters') +
  
  labs(x = "UMAP1", y = "UMAP2") +
  
  theme(axis.text.y = element_blank(),
        
        text = element_text(size = 20),
        
        axis.ticks.y = element_blank(),
        
        axis.text.x = element_blank(),
        
        axis.ticks.x = element_blank())+
  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

## Fig. 6B
bubble.df=t(mega_seurat@assays$RNA@scale.data)

bubble.df=bubble.df[,marker]

bubble.df <- as.data.frame(bubble.df)

bubble.df$CB=rownames(bubble.df)

mega_seurat_metadata <- mega_seurat@meta.data

mega_seurat_metadata$CB <- rownames(mega_seurat_metadata)

bubble.df=merge(bubble.df,mega_seurat_metadata[,c("CB",'seurat_clusters')],by = "CB")

celltype_v=c()

gene_v=c()

mean_v=c()

ratio_v=c()

for (i in unique(bubble.df$seurat_clusters)) {
  
  bubble.df_small=bubble.df%>%filter(seurat_clusters==i)
  
  for (j in marker) {
    
    exp_mean=mean(bubble.df_small[,j])
    
    exp_ratio=sum(bubble.df_small[,j] > min(bubble.df_small[,j])) / length(bubble.df_small[,j])
    
    celltype_v=append(celltype_v,i)
    
    gene_v=append(gene_v,j)
    
    mean_v=append(mean_v,exp_mean)
    
    ratio_v=append(ratio_v,exp_ratio)
    
  }
}
plotdf=data.frame(
  
  celltype=celltype_v,
  
  gene=gene_v,
  
  exp=mean_v,
  
  ratio=ratio_v
)
plotdf$celltype=factor(plotdf$celltype,levels = sort(unique(plotdf$celltype)))

plotdf$gene=factor(plotdf$gene,levels = rev(as.character(marker)))

plotdf$exp=ifelse(plotdf$exp>3,3,plotdf$exp)

plotdf <- plotdf[-c(which(plotdf$ratio==0)),]


plotdf%>%ggplot(aes(x=celltype,y=gene,size=ratio,color=exp))+geom_point()+
  
  scale_x_discrete("")+scale_y_discrete("")+
  
  scale_color_gradientn(colours = c('#8DB59C','#A9C0A9','#E9E6C5','#DFC3A5','#D68B7D','#D47679'))+
  
  scale_size_continuous(limits = c(0,1))+theme_bw()+scale_fill_continuous(breaks=c(0,0.5,1))+
  
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 55),axis.text.x = element_text(size = 10,colour = "black"),
    
    axis.text.y = element_text(size = 10,colour = "black")
  )

## Fig. 6E

ggplot(mega_gsva_score,aes(x=seurat_cluster, y=Inositol_phosphate_metabolism,color=seurat_cluster,fill=seurat_cluster,alpha=0.6))+
  
  geom_boxplot(color = c(cell_type_cols[1:7]),size=1,width = 0.5,outlier.colour = NA) + theme_bw()  + 
  
  theme(panel.grid = element_blank(),
        
        panel.background = element_blank(),
        
        axis.line = element_line(),
        
        axis.text.x = element_text(size=12,color="black", face= "bold",angle=60,hjust = 1),
        
        axis.text.y = element_text(size=12, color="black", face= "bold"),
        
        plot.title = element_text(size=12,color="black", face= "bold",hjust=0))+
  
  #geom_jitter(width = 0.2,show.legend =F,size=1.2) +
  scale_fill_manual(values=c(cell_type_cols[1:7]))+
  
  labs(x="",y="GSVA score",title="Inositol_phosphate_metabolism")+#ylim(0.1,0.3)+
  
  theme(strip.text.x = element_text(size = 12, colour = "black"))+
  
  theme(strip.background.x = element_rect(fill = "lightgrey")) 


mydata=FetchData(mega_seurat,vars = c(
  'UMAP_1','UMAP_2'))

as.numeric(mega_gsva_score[,'Inositol_phosphate_metabolism']) -> ip


mydata$selg <- ip

colnames(mydata) <- gsub(' ','_',colnames(mydata))


mydata <- mydata[order(mydata$selg),]   

p <- ggplot(mydata,aes(x=UMAP_1,y=UMAP_2))+
  
  geom_point(data = mydata,aes(x=UMAP_1,y=UMAP_2,
                               
                               color=selg),size=0.1)+
  
  scale_colour_gradientn('selg',limits=c(0,0.2),colors = c('#EAEBEB','yellow','red','darkred'))+
  
  theme_bw()+theme(panel.grid.major=element_line(colour=NA),
                   
                   panel.background = element_rect(fill = "transparent",colour = NA),
                   
                   plot.background = element_rect(fill = "transparent",colour = NA),
                   
                   panel.grid.minor = element_blank())
p


### Fig.6F
data <- as(as.matrix(mega_seurat@assays$RNA@counts),"sparseMatrix")
p_data <- mega_seurat@meta.data

f_data <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

pd <- new("AnnotatedDataFrame", data = p_data)

fd <- new("AnnotatedDataFrame", data = f_data)

cds <- newCellDataSet(data,
                      
                      phenoData = pd,
                      
                      featureData = fd,
                      
                      lowerDetectionLimit = 0.5,
                      
                      expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)

cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)

expressed_genes <- VariableFeatures(mega_seurat)

cds <- setOrderingFilter(cds, expressed_genes)

plot_ordering_genes(cds)

diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~seurat_clusters", cores = 1) 

deg <- subset(diff, qval < 0.01)

deg <- deg[order(deg$qval,decreasing = F),]

head(deg)
ordergene <- rownames(deg)

cds <- setOrderingFilter(cds,ordergene)

plot_ordering_genes(cds) 

ordergene <- row.names(deg)[order(deg$pval)][1:400] 

cds <- reduceDimension(cds, max_components = 2, method = "DDRTree")

cds <- orderCells(cds)

plot_cell_trajectory(cds, color_by = "sct_cluster_new",size=1, show_backbone = TRUE)+scale_color_manual(values = cell_type_cols)


### Fig. 6G
Time_diff <- differentialGeneTest(cds, cores = 1,
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")

Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)] 

Time_diff_genes <- Time_diff[c(which(Time_diff$pval < 0.05)),]

Time_diff_genes <- Time_diff_genes[c(which(Time_diff_genes$qval < 0.05)),]

Time_diff_genes <- Time_diff_genes[!grepl('ncRNA',Time_diff_genes$gene_short_name),] 

Time_diff_genes <- Time_diff_genes[grepl('-',Time_diff_genes$gene_short_name),]

Time_diff_genes <- Time_diff_genes[order(Time_diff_genes$pval),]

Time_genes <- Time_diff_genes %>% pull(gene_short_name) %>% as.character()
#Time_genes <- Time_diff_genes$gene_short_name[1:100]

Time_genes <- Time_diff_genes$gene_short_name
dev.off()

p <- plot_pseudotime_heatmap(cds[Time_genes,],num_clusters=3, show_rownames = F, return_heatmap = T,hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))

clusters <- cutree(p$tree_row,k = 3)

clusters <- as.data.frame(clusters)

clusters[,1] <- as.character(clusters[,1])

colnames(clusters) <- "gene_clusters"

gmtFile_mega <- "gmt_mega.gmt"

geneSets_mega_total <- getGmt(gmtFile_mega) #signature read

id <- c()

for (i in 1:length(geneSets)){
  
  c(geneSets[[i]]@geneIds) -> gene_c
  
  gene_c <- gsub('_','-',gene_c)
  
  length(gene_c) -> ip
  
  gene_c <- as.data.frame(gene_c)
  
  gene_c$ID <- geneSets[[i]]@shortDescription
  
  gene_c$Description <- geneSets[[i]]@setName
  
  id <- rbind(id,gene_c)
}
id <- as.data.frame(id)

colnames(id) <- c('gene_id','ID','Description')

table(clusters$cluster)

time_gene_cluster1 <- rownames(clusters)[which(clusters$cluster=="1")]

time_gene_cluster2 <- rownames(clusters)[which(clusters$cluster=="2")]

time_gene_cluster3 <- rownames(clusters)[which(clusters$cluster=="3")]

go_rich_1 <- enricher(gene = rownames(clusters),  
                    
                    TERM2GENE = id[c('ID', 'gene_id')], 
                    
                    TERM2NAME = id[c('ID', 'Description')], 
                    
                    pAdjustMethod = 'BH', 
                    
                    pvalueCutoff = 1,  
                    
                    qvalueCutoff = 1)  


go_rich_1 <- enricher(gene = time_gene_cluster1,
                    
                    TERM2GENE = id[c('ID', 'gene_id')],  
                    
                    TERM2NAME = id[c('ID', 'Description')],pvalueCutoff = 1,minGSSize = 1,
                    
                    qvalueCutoff = 1) 


go_rich_2 <- enricher(gene = rownames(clusters),  
                      
                      TERM2GENE = id[c('ID', 'gene_id')], 
                      
                      TERM2NAME = id[c('ID', 'Description')], 
                      
                      pAdjustMethod = 'BH', 
                      
                      pvalueCutoff = 1,  
                      
                      qvalueCutoff = 1)  


go_rich_2 <- enricher(gene = time_gene_cluster2,
                      
                      TERM2GENE = id[c('ID', 'gene_id')],  
                      
                      TERM2NAME = id[c('ID', 'Description')],pvalueCutoff = 1,minGSSize = 1,
                      
                      qvalueCutoff = 1) 

go_rich_3 <- enricher(gene = rownames(clusters),  
                      
                      TERM2GENE = id[c('ID', 'gene_id')], 
                      
                      TERM2NAME = id[c('ID', 'Description')], 
                      
                      pAdjustMethod = 'BH', 
                      
                      pvalueCutoff = 1,  
                      
                      qvalueCutoff = 1)  


go_rich_3 <- enricher(gene = time_gene_cluster3,
                      
                      TERM2GENE = id[c('ID', 'gene_id')],  
                      
                      TERM2NAME = id[c('ID', 'Description')],pvalueCutoff = 1,minGSSize = 1,
                      
                      qvalueCutoff = 1) 

ggplot(go_rich_1,aes(y=Description,x=Count))+
  
  geom_bar(stat = "identity",width=0.7,fill = "#4DBBD5")+

  labs(title = "",
       
       x = "Gene numbers", 
       
       y = "Pathways")+
  
  theme_bw()+
  
  theme(axis.title.x = element_text(face = "bold",size = 8),
        
        axis.title.y = element_text(face = "bold",size = 8),
        
        legend.title = element_text(face = "bold",size = 8),
        
        panel.grid.major=element_line(colour=NA),
        
        panel.grid.minor = element_blank())


ggplot(go_rich_2,aes(y=Description,x=Count,fill=pvalue))+
  
  geom_bar(stat = "identity",position = "dodge")+

  labs(title = "",
       x = "Gene numbers", 
       
       y = "Pathways")+
  
  theme_bw()+
  
  scale_fill_gradient(high = '#91D1C2',low = '#F3DA7D')+
  
  theme(axis.title.x = element_text(face = "bold",size = 8),
        
        axis.title.y = element_text(face = "bold",size = 8),
        
        legend.title = element_text(face = "bold",size = 8),
        
        panel.grid.major=element_line(colour=NA),
        
        panel.grid.minor = element_blank())

ggplot(go_rich_3,aes(y=Description,x=Count,fill=pvalue))+
  
  geom_bar(stat = "identity",position = "dodge")+

  labs(title = "",
       
       x = "Gene numbers", 
       
       y = "Pathways")+
  theme_bw()+
  scale_fill_gradient(high = '#71A36D',low = '#F3DA7D')+
  
  theme(axis.title.x = element_text(face = "bold",size = 8),
        
        axis.title.y = element_text(face = "bold",size = 8),
        
        legend.title = element_text(face = "bold",size = 8),
        
        panel.grid.major=element_line(colour=NA),
        
        panel.grid.minor = element_blank())

###
