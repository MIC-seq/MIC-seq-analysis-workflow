### Figure 2 ##########
library(ggplot2)
library(Seurat)
#### Fig.2 A
seurat_object <- readRDS('combine_dataset.rds')

DimPlot(seurat_object, reduction = "umap",label.size = 6,label = F,
        
        group.by = 'species_info',cols = cell_type_cols) +
  
  #NoLegend() +
  labs(x = "UMAP1", y = "UMAP2") +
  
  theme(axis.text.y = element_blank(),
        
        text = element_text(size = 20),
        
        axis.ticks.y = element_blank(),
        
        axis.text.x = element_blank(),
        
        #axis.title = element_text(size = 20),
        
        axis.ticks.x = element_blank())+
  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

### Fig. 2B
tree <- parseMetaphlanTSV(cbx, node.size.offset=1, node.size.scale=0.8)

p <- tree.backbone(tree, size=0.5)
p
tree@data -> tree_metadata

tree@data$Phylum<- sku[1:81]

ggtree(tree, aes(Phylum))+  theme(legend.position = "right")+geom_tiplab(size=3)

gra <- ggtree(tree, aes(color=Phylum), layout="circular", size=0.1)  +
  
  geom_tiplab(aes(label=label, col=Phylum), hjust=-0.5, align=TRUE, linesize=0.2)+
  
  theme(legend.title=element_text(face="bold"), legend.position="bottom", legend.box="horizontal", legend.text=element_text(size=rel(0.5)))
gra

### Fig. 2C
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
plotdf <- dotplot_analysis(sample = seurat_object,marker = choose_marker_gene,colnames_show = choose_marker_species$KEGGnumber,sample_metadata = seurat_object@meta.data,select_info='species_info')

plotdf%>%ggplot(aes(x=celltype,y=gene,size=ratio,color=exp))+geom_point()+
  
  scale_x_discrete("")+scale_y_discrete("")+
  
  scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11])))+
  
  scale_size_continuous(limits = c(0,1))+theme_bw()+scale_fill_continuous(breaks=c(0,0.5,1))+
  
  
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45),axis.text.x = element_text(size = 10,colour = "black"),
    axis.text.y = element_text(size = 10,colour = "black")
  )

### Fig. 2D

pdata <- rbind(species_donor,specis_anno)

nodes <- data.frame(name = c(as.character(pdata$source), 
                             as.character(pdata$target)) %>% unique())

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

### Fig. 2E
seurat_gene_umap = FetchData(seurat_object,vars = c(
  'UMAP_1','UMAP_2'))

as.matrix(seurat_object@assays$RNA@counts) -> gene_matrix
select_genename <- "aprA"

as.numeric(gene_matrix$select_genename) -> ip

seurat_gene_umap$select_gene <- ip

ggplot(seurat_gene_umap,aes(x=UMAP_1,y=UMAP_2))+
  
  geom_point(data = seurat_gene_umap,aes(x=UMAP_1,y=UMAP_2,
                                         
                               color=select_gene),size=0.1)+
  
  scale_colour_gradientn('select_gene',colors = c('#EAEBEB','yellow','red','darkred'))+
  
  theme_bw()+
  
  theme(panel.grid.major=element_line(colour=NA),
        
                   panel.background = element_rect(fill = "transparent",colour = NA),
        
                   plot.background = element_rect(fill = "transparent",colour = NA),
        
                   panel.grid.minor = element_blank())


### Fig. 2F
seurat_gene_umap = FetchData(seurat_object,vars = c(
  'UMAP_1','UMAP_2'))

as.matrix(seurat_object@assays$RNA@counts) -> gene_matrix
select_genename <- "hag"

as.numeric(gene_matrix$select_genename) -> ip

seurat_gene_umap$select_gene <- ip

ggplot(seurat_gene_umap,aes(x=UMAP_1,y=UMAP_2))+
  
  geom_point(data = seurat_gene_umap,aes(x=UMAP_1,y=UMAP_2,
                                         
                                         color=select_gene),size=0.1)+
  
  scale_colour_gradientn('select_gene',colors = c('#EAEBEB','yellow','red','darkred'))+
  
  theme_bw()+
  
  theme(panel.grid.major=element_line(colour=NA),
        
        panel.background = element_rect(fill = "transparent",colour = NA),
        
        plot.background = element_rect(fill = "transparent",colour = NA),
        
        panel.grid.minor = element_blank())





### Fig. 2G
seurat_gene_umap = FetchData(seurat_object,vars = c(
  'UMAP_1','UMAP_2'))

as.matrix(seurat_object@assays$RNA@counts) -> gene_matrix
select_genename <- "sodB"

as.numeric(gene_matrix$select_genename) -> ip

seurat_gene_umap$select_gene <- ip

ggplot(seurat_gene_umap,aes(x=UMAP_1,y=UMAP_2))+
  
  geom_point(data = seurat_gene_umap,aes(x=UMAP_1,y=UMAP_2,
                                         
                                         color=select_gene),size=0.1)+
  
  scale_colour_gradientn('select_gene',colors = c('#EAEBEB','yellow','red','darkred'))+
  
  theme_bw()+
  
  theme(panel.grid.major=element_line(colour=NA),
        
        panel.background = element_rect(fill = "transparent",colour = NA),
        
        plot.background = element_rect(fill = "transparent",colour = NA),
        
        panel.grid.minor = element_blank())

