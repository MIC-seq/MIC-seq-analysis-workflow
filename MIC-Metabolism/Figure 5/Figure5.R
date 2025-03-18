########Figure 5

## Fig. 5A

DimPlot(seurat_etp, reduction = "umap",label.size = 6,label = F,cols = c("yellow","organge",'grey'),
        
        group.by = 'tp') +
  
  #NoLegend() +
  labs(x = "UMAP1", y = "UMAP2") +
  
  theme(axis.text.y = element_blank(),
        
        text = element_text(size = 20),
        
        axis.ticks.y = element_blank(),
        
        axis.text.x = element_blank(),
        
        #axis.title = element_text(size = 20),
        axis.ticks.x = element_blank())+
  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))


DimPlot(seurat_etp, reduction = "umap",label.size = 6,label = F,cols = cell_type_cols,
        
        group.by = 'species_info') +
  #NoLegend() +
  labs(x = "UMAP1", y = "UMAP2") +
  
  theme(axis.text.y = element_blank(),
        
        text = element_text(size = 20),
        
        axis.ticks.y = element_blank(),
        
        axis.text.x = element_blank(),
        
        #axis.title = element_text(size = 20),
        axis.ticks.x = element_blank())+
  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

DimPlot(seurat_etf, reduction = "umap",label.size = 6,label = F,cols = c("yellow","organge",'grey'),
        
        group.by = 'tp') +
  
  #NoLegend() +
  labs(x = "UMAP1", y = "UMAP2") +
  
  theme(axis.text.y = element_blank(),
        
        text = element_text(size = 20),
        
        axis.ticks.y = element_blank(),
        
        axis.text.x = element_blank(),
        
        #axis.title = element_text(size = 20),
        axis.ticks.x = element_blank())+
  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))


DimPlot(seurat_etf, reduction = "umap",label.size = 6,label = F,cols = cell_type_cols,
        
        group.by = 'species_info') +
  #NoLegend() +
  labs(x = "UMAP1", y = "UMAP2") +
  
  theme(axis.text.y = element_blank(),
        
        text = element_text(size = 20),
        
        axis.ticks.y = element_blank(),
        
        axis.text.x = element_blank(),
        
        #axis.title = element_text(size = 20),
        axis.ticks.x = element_blank())+
  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

DimPlot(seurat_etb, reduction = "umap",label.size = 6,label = F,cols = c("yellow","organge",'grey'),
        
        group.by = 'tp') +
  
  #NoLegend() +
  labs(x = "UMAP1", y = "UMAP2") +
  
  theme(axis.text.y = element_blank(),
        
        text = element_text(size = 20),
        
        axis.ticks.y = element_blank(),
        
        axis.text.x = element_blank(),
        
        #axis.title = element_text(size = 20),
        axis.ticks.x = element_blank())+
  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))


DimPlot(seurat_etb, reduction = "umap",label.size = 6,label = F,cols = cell_type_cols,
        
        group.by = 'species_info') +
  #NoLegend() +
  labs(x = "UMAP1", y = "UMAP2") +
  
  theme(axis.text.y = element_blank(),
        
        text = element_text(size = 20),
        
        axis.ticks.y = element_blank(),
        
        axis.text.x = element_blank(),
        
        #axis.title = element_text(size = 20),
        axis.ticks.x = element_blank())+
  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

### Fig. 5c

pdata <- rbind(layer_12,layer_23)
nodes <- data.frame(name = c(as.character(pdata$source), 
                             as.character(pdata$target)) %>% unique())
# 把source、target转换为数字
pdata$IDsource = match(pdata$source, nodes$name)-1

pdata$IDtarget = match(pdata$target, nodes$name)-1

head(pdata)

library(networkD3)
library(webshot)
p1 <- sankeyNetwork(Links = pdata, 
                    
                    Nodes = nodes, 
                    
                    Source = "IDsource", 
                    
                    
                    Target = "IDtarget", 
                    
                    Value = "weight",
                    
                    NodeID = "name",
                    
                    LinkGroup = 'source', 
                    
                    sinksRight = FALSE, 
                    
                    nodeWidth = 5, 
                    
                    fontSize = 18, 
                    
                    nodePadding = 4) 

p1


### Fig. 5D

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
 
  datac <- plyr::rename(datac, c("mean" = measurevar))
  

  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

tgc <- summarySE(d1tp123_each_barcode_zscore_function, measurevar="len", groupvars=c("supp","dose"))

tgc <- as.data.frame(tgc)

tgc$dose <- as.numeric(tgc$dose)

ggplot(tgc, aes(x=dose, y=len, colour=supp)) + 
  
  geom_errorbar(aes(ymin=len-se, ymax=len+se), width=.1) + geom_point() +
  
  geom_line(group = 1) +
  
  theme_bw() + scale_color_manual(values = cell_type_cols) 

ggplot(data=path_ca, aes(x=variable, y=value,color = species)) + geom_point() + geom_line(group = 1) + 
  
  geom_errorbar(aes(ymax = value + sd, ymin = value -  sd),width = 0.15) +
  
  facet_grid(path~species) +theme_bw() + scale_color_manual(values = cell_type_cols) 

tgc <- summarySE(d2tp123_each_barcode_zscore_function, measurevar="len", groupvars=c("supp","dose"))

tgc <- as.data.frame(tgc)

tgc$dose <- as.numeric(tgc$dose)

ggplot(tgc, aes(x=dose, y=len, colour=supp)) + 
  
  geom_errorbar(aes(ymin=len-se, ymax=len+se), width=.1) + geom_point() +
  
  geom_line(group = 1) +
  
  theme_bw() + scale_color_manual(values = cell_type_cols) 

ggplot(data=path_ca, aes(x=variable, y=value,color = species)) + geom_point() + geom_line(group = 1) + 
  
  geom_errorbar(aes(ymax = value + sd, ymin = value -  sd),width = 0.15) +
  
  facet_grid(path~species) +theme_bw() + scale_color_manual(values = cell_type_cols) 

tgc <- summarySE(d3tp123_each_barcode_zscore_function, measurevar="len", groupvars=c("supp","dose"))

tgc <- as.data.frame(tgc)

tgc$dose <- as.numeric(tgc$dose)

ggplot(tgc, aes(x=dose, y=len, colour=supp)) + 
  
  geom_errorbar(aes(ymin=len-se, ymax=len+se), width=.1) + geom_point() +
  
  geom_line(group = 1) +
  
  theme_bw() + scale_color_manual(values = cell_type_cols) 

ggplot(data=path_ca, aes(x=variable, y=value,color = species)) + geom_point() + geom_line(group = 1) + 
  
  geom_errorbar(aes(ymax = value + sd, ymin = value -  sd),width = 0.15) +
  
  facet_grid(path~species) +theme_bw() + scale_color_manual(values = cell_type_cols) 

### Fig.5F
pc.num <- 1:30

clustering_analysis <- function(sample=NA,resolution_num=NA){
  
  sample <- NormalizeData(sample, normalization.method = "LogNormalize", scale.factor = 10000)
  
  #sample <- SCTransform(sample)
  sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 3000)
  
  sample <- ScaleData(sample,features = VariableFeatures(sample))
  
  sample <- RunPCA(sample,features = VariableFeatures(sample))
  
  #sample <- sample %>% RunHarmony("orig.ident", plot_convergence = TRUE, lambda = 0.5)
  sample <- FindNeighbors(sample,dims = pc.num, reduction = "pca")
  
  sample <- FindClusters(sample, resolution = resolution_num)
  
  sample <- RunUMAP(sample, dims = pc.num, reduction = "pca")
  return(sample)
}

dotplot_analysis <- function(sample = NA,marker=NA,colnames_show=NA,sample_metadata = NA,select_info=NA){
  
  bubble.df=t(sample@assays$RNA@scale.data)
  
  bubble.df=bubble.df[,marker]
  
  bubble.df <- as.data.frame(bubble.df)
  
  colnames(bubble.df) <- colnames_show
  
  bubble.df$CB=rownames(bubble.df)
  
  bubble.df=merge(bubble.df,sample_metadata[,c("CB",select_info)],by = "CB")
  
  celltype_v=c()
  
  gene_v=c()
  
  mean_v=c()
  
  ratio_v=c()
  
  for (i in unique(bubble.df$species_info)) {
    
    bubble.df_small=bubble.df%>%filter(species_info==i)
    
    for (j in colnames_show) {
      
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
  
  plotdf$gene=factor(plotdf$gene,levels = rev(as.character(colnames_show)))
  
  plotdf$exp=ifelse(plotdf$exp>3,3,plotdf$exp)
  
  plotdf <- plotdf[-c(which(plotdf$ratio==0)),]
  
  return(plotdf)
}

pcopri_seurat <- clustering_analysis(pcopri_seurat,resolution_num = 0.3)

pcopri_seurat@meta.data -> pcopri_metadata

DimPlot(pcopri_seurat, reduction = "umap",label.size = 6,label = F,cols = cell_type_cols,
        
        group.by = 'seurat_clusters') +
  
  #NoLegend() +
  
  labs(x = "UMAP1", y = "UMAP2") +
  
  theme(axis.text.y = element_blank(),
        
        text = element_text(size = 20),
        
        axis.ticks.y = element_blank(),
        
        axis.text.x = element_blank(),
        
        #axis.title = element_text(size = 20),
        axis.ticks.x = element_blank())+
  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

pcopri_seurat@meta.data$anno1 <- pcopri_seurat@meta.data$orig.ident

tab.1=table(pcopri_seurat$seurat_clusters,pcopri_seurat$anno1) 

balloonplot(tab.1)

Idents(pcopri_seurat) <- pcopri_seurat$seurat_clusters

pcopri_marker <- FindAllMarkers(pcopri_seurat,logfc.threshold = 0.01,min.pct = 0.01)

plotdf <- dotplot_analysis(sample = pcopri_seurat,marker = choose_marker_gene,colnames_show = choose_marker_species$KEGGnumber,sample_metadata = pcopri_seurat@meta.data,select_info='seurat_clusters')

plotdf%>%ggplot(aes(x=celltype,y=gene,size=ratio,color=exp))+geom_point()+
  
  scale_x_discrete("")+scale_y_discrete("")+
  
  scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11])))+
  
  scale_size_continuous(limits = c(0,1))+theme_bw()+scale_fill_continuous(breaks=c(0,0.5,1))+
  
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45),axis.text.x = element_text(size = 10,colour = "black"),
    axis.text.y = element_text(size = 10,colour = "black")
  )

## Fig. 5G

b_seurat <- clustering_analysis(b_seurat,resolution_num = 0.3)

b_seurat@meta.data -> b_metadata

DimPlot(b_seurat, reduction = "umap",label.size = 6,label = F,cols = cell_type_cols,
        
        group.by = 'seurat_clusters') +
  
  #NoLegend() +
  
  labs(x = "UMAP1", y = "UMAP2") +
  
  theme(axis.text.y = element_blank(),
        
        text = element_text(size = 20),
        
        axis.ticks.y = element_blank(),
        
        axis.text.x = element_blank(),
        
        #axis.title = element_text(size = 20),
        axis.ticks.x = element_blank())+
  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

b_seurat@meta.data$anno1 <- b_seurat@meta.data$orig.ident

tab.1=table(b_seurat$seurat_clusters,b_seurat$anno1) 

balloonplot(tab.1)

Idents(b_seurat) <- b_seurat$seurat_clusters

b_marker <- FindAllMarkers(b_seurat,logfc.threshold = 0.01,min.pct = 0.01)

plotdf <- dotplot_analysis(sample = b_seurat,marker = choose_marker_gene,colnames_show = choose_marker_species$KEGGnumber,sample_metadata = b_seurat@meta.data,select_info='seurat_clusters')

plotdf%>%ggplot(aes(x=celltype,y=gene,size=ratio,color=exp))+geom_point()+
  
  scale_x_discrete("")+scale_y_discrete("")+
  
  scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11])))+
  
  scale_size_continuous(limits = c(0,1))+theme_bw()+scale_fill_continuous(breaks=c(0,0.5,1))+
  
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45),axis.text.x = element_text(size = 10,colour = "black"),
    axis.text.y = element_text(size = 10,colour = "black")
  )






