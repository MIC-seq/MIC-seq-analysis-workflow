#!/usr/bin/R

library(dplyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(Seurat)
library(tidyverse) 


###Single-cell transcriptional analysis for host-phage relationships
cluster_pipeline <- function(input_data=NA, bacteria_info = NA, min_cells=30, min_features=15,
 sp_rate_threshold = 0.03, pc.num=1:30, resolution_num = 0.5, if_nr = FALSE) {
  ##Seurat对象创建
  sample <- CreateSeuratObject(counts = input_data, project = "A197",min.cells = min_cells, min.features = min_features)
  ##物种信息导入
  BC <- colnames(sample)
  BC<-as.data.frame(BC)
  colnames(bacteria_info) <- c("BC","index")
  BC_info <- merge(BC, bacteria_info,all.x=TRUE)
  sample@meta.data$genus_info <- BC_info[,2]
  genus_num_data = table(sample@meta.data$genus_info)

  ##去除极低丰度菌(<3%)
  num_threshold = dim(sample)[2] * sp_rate_threshold
  for(i in 1:length(genus_num_data)){
    if(genus_num_data[[i]]<num_threshold){
      index = names(genus_num_data[i])
      sample <- subset(sample,subset = genus_info != index)
    }
  }

  ##样本基因UMI小提琴图
  p1<-VlnPlot(sample, features = "nFeature_RNA", col="green",pt.size = 0)+
    geom_boxplot(width=.3, col="black", fill="white")+
    NoLegend()

  p2 <- VlnPlot(sample, features = "nCount_RNA", col="blue",pt.size = 0)+
    geom_boxplot(width=.3, col="black", fill="white")+
    NoLegend()

  ##data UMI gene count
  if(if_nr == TRUE){
    filename1 = paste("features_nr_",min_cells,"_",min_features,".pdf", sep = "")
  }else{
    filename1 = paste("features_",min_cells,"_",min_features,".pdf", sep = "")
  }
  ggsave(filename = filename1,
        plot = p1|p2,
        width = 10,
        height = 7,
        units = "in",
        dpi=500)

  ##降维聚类分析
  sample <- NormalizeData(sample, normalization.method = "LogNormalize", scale.factor = 10000)
  sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000)
  sample <- ScaleData(sample,features = VariableFeatures(sample))
  sample <- RunPCA(sample,features = VariableFeatures(sample))
  sample <- FindNeighbors(sample,dims = pc.num)
  sample <- FindClusters(sample, resolution = resolution_num)
  sample <- RunUMAP(sample, dims = pc.num)

  ##提取UMAP聚类数据信息
  umap = sample@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(seurat_clusters = sample@meta.data$seurat_clusters) %>% 
  cbind(genus_info = sample@meta.data$genus_info)

  seurat_clusters_med <- umap %>%
  group_by(seurat_clusters) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

  genus_info_med <- umap %>%
  group_by(genus_info) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

  cell_type_cols <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#F08080","#1E90FF",
                      "#808000","#FF00FF","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
                      "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
                      "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

  plot1 = DimPlot(sample, reduction = "umap", label.size = 6 , label = T, pt.size = 1.5, cols= cell_type_cols, group.by = 'seurat_clusters') +
                  NoLegend() +
                  labs(x = "UMAP1", y = "UMAP2") +
                  theme(axis.text.y = element_blank(),
                  text = element_text(size = 20,face="bold"),
                  axis.ticks.y = element_blank(),
                  axis.text.x = element_blank(),
                  axis.title = element_text(size = 20),
                  axis.ticks.x = element_blank()) +
                  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) 
                  


  plot2 <- ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = genus_info)) + 
            geom_point(size = 1.5 , alpha =1 ) + 
            scale_color_manual(values = cell_type_cols) +
            ggtitle("genus_info") +
            labs(x = "UMAP_1", y = "UMAP_2") +
            theme(plot.title = element_text(hjust = 0.5), #标题居中
            panel.grid.major = element_blank(), #主网格线
            panel.grid.minor = element_blank(), #次网格线
            panel.border = element_blank(), #边框
            #axis.title = element_blank(),  #轴标题
            axis.text = element_blank(), # 坐标刻度
            text = element_text(size = 20),
            axis.ticks = element_blank(),
            panel.background = element_rect(fill = 'white'), #背景色
            plot.background=element_rect(fill="white")) +
            theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) +
            geom_label_repel(aes(label=genus_info), fontface="bold",data = genus_info_med,
            point.padding=unit(0.5, "lines")) +
            NoLegend()

  ##物种信息与UMAP聚类图整合展示
  if(if_nr == TRUE){
    filename2 = paste("test_nr_",min_cells,"_",min_features,"_",resolution_num,".pdf", sep = "")
  }else{
    filename2 = paste("test_",min_cells,"_",min_features,"_",resolution_num,".pdf", sep = "")
  }
  ggsave(filename = filename2,
        plot = plot1+plot2,
        width = 20,
        height = 10,
        units = "in",
        dpi=500)

  return(sample)
}
sample <- cluster_pipeline(input_data = VC_nonrna_data_A197, bacteria_info = genus_info_A197, min_cells = 15, min_features = 15)


###Filter contaminated microbes
sample_fig <- subset(sample, seurat_clusters != 0 & nFeature_RNA <100)
sample_fig <- subset(sample_fig, genus_info != "CAG-81" | seurat_clusters == 9 | seurat_clusters == 13)
sample_fig <- subset(sample_fig, genus_info != "Clostridium_Q " | seurat_clusters == 3 | seurat_clusters == 14)
sample_fig <- subset(sample_fig, genus_info != "Dorea_A" | seurat_clusters == 8)
sample_fig <- subset(sample_fig, genus_info != "Faecalibacterium" | seurat_clusters == 4)
sample_fig <- subset(sample_fig, genus_info != "Fusicatenibacter" | seurat_clusters == 6)
sample_fig <- subset(sample_fig, genus_info != "Lachnospira" | seurat_clusters == 11 | seurat_clusters == 15)
sample_fig <- subset(sample_fig, genus_info != "Phascolarctobacterium_A" | seurat_clusters == 2)
sample_fig <- subset(sample_fig, genus_info != "Prevotella" | seurat_clusters == 1| seurat_clusters == 10 )
sample_fig <- subset(sample_fig, genus_info != "Roseburia" | seurat_clusters == 5| seurat_clusters == 7 | seurat_clusters == 12)
redraw <- function(sample = sample_fig, pc.num = 1:30, resolution_num = 0.5){
  sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000)
  sample <- ScaleData(sample,features = VariableFeatures(sample))
  sample <- RunPCA(sample,features = VariableFeatures(sample))
  sample <- FindNeighbors(sample,dims = pc.num)
  sample <- FindClusters(sample, resolution = resolution_num)
  sample <- RunUMAP(sample, dims = pc.num)

  umap = sample@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(seurat_clusters = sample@meta.data$seurat_clusters) %>% 
  cbind(genus_info = sample@meta.data$genus_info)

  seurat_clusters_med <- umap %>%
  group_by(seurat_clusters) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

  genus_info_med <- umap %>%
  group_by(genus_info) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

  cell_type_cols <- c(brewer.pal(9, "Set1"), "#FF34B3","#BC8F8F","#20B2AA","#00F5FF",
                    "#FFA500","#ADFF2F","#FF6A6A","#7FFFD4", "#AB82FF","#90EE90",
                    "#00CD00","#008B8B","#6495ED","#FFC1C1","#CD5C5C","#8B008B",
                    "#FF3030", "#7CFC00","#000000","#708090")

  plot1 = DimPlot(sample, reduction = "umap", label.size = 9 , label = T, pt.size = 1, cols= cell_type_cols, group.by = 'seurat_clusters') +
                  NoLegend() +
                  labs(x = "UMAP_1", y = "UMAP_2") +
                  theme(axis.text.y = element_blank(),
                  text = element_text(size = 25),
                  axis.ticks.y = element_blank(),
                  axis.text.x = element_blank(),
                  axis.title = element_text(size = 25),
                  axis.ticks.x = element_blank()) +
                  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) 
                  


  plot2 <- ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = genus_info)) + 
            geom_point(size = 1 , alpha =1 ) + 
            scale_color_manual(values = cell_type_cols) +
            ggtitle("genus_info") +
            labs(x = "UMAP_1", y = "UMAP_2") +
            theme(plot.title = element_text(hjust = 0.5), #标题居中
            panel.grid.major = element_blank(), #主网格线
            panel.grid.minor = element_blank(), #次网格线
            panel.border = element_blank(), #边框
            #axis.title = element_blank(),  #轴标题
            axis.text = element_blank(), # 坐标刻度
            text = element_text(size = 25),
            axis.ticks = element_blank(),
            panel.background = element_rect(fill = 'white'), #背景色
            plot.background=element_rect(fill="white")) +
            theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) +
            geom_label_repel(aes(label=genus_info),data = genus_info_med, direction = "both",
            size=8,nudge_x = -2 , nudge_y = 2.5, 
            point.padding=unit(0.5, "lines")) +
            NoLegend()

  ggsave(filename = "aaa.pdf",
        plot = plot1+plot2,
        width = 15,
        height = 7,
        units = "in",
        dpi=500)

  return(sample)
}
sample_fig <- redraw(sample_fig)

new.cluster.ids <- c("Prevotella","Phascolarcto","Clostridium","Faecalibacterium",
                    "Roseburia","Fusicatenibacter","Roseburia","Dorea","CAG-81","Prevotella",
                    "Lachnospira","Roseburia","CAG-81","Clostridium","Lachnospira")
names(new.cluster.ids) <- levels(sample_fig)
sample_fig <- RenameIdents(sample_fig, new.cluster.ids)
Idents(sample_fig) <- factor(Idents(sample_fig),levels = c("Prevotella", "Phascolarcto",
                      "Clostridium","Dorea","Roseburia","Lachnospira","Faecalibacterium",
                      "CAG-81","Fusicatenibacter"))


###phage numbers count for each genus
plot <- VlnPlot(sample_fig, features = "nFeature_RNA", col=cell_type_cols,pt.size = 0)+ 
        NoLegend()+
        geom_boxplot(width=.1, col="black", lwd = 1.2, fatten = 1.5, outlier.shape = NA) + 
        theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 30),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=30)) + 
        theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

ggsave(filename = "nFeature_virus_cluster.pdf",
      plot = plot,
      width = 15,
      height = 7,
      units = "in",
      dpi=500)


###major transcriptional host-phage relationship analysis
sample_markers <- FindAllMarkers(sample_fig, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1) 

sample_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

sample_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5

Prevotella_markers = FindMarkers(sample_fig,ident.1 = "Prevotella", min.pct = 0.1)
Phascolarctobacterium_markers = FindMarkers(sample_fig, ident.1 = "Phascolarctobacterium_A", min.pct = 0.1)
Dorea_markers = FindMarkers(sample_fig, ident.1 = "Dorea_A", min.pct = 0.1)
Clostridium_markers = FindMarkers(sample_fig, ident.1 = "Clostridium_Q", min.pct = 0.1)
Roseburia_markers = FindMarkers(sample_fig, ident.1 = "Roseburia", min.pct = 0.1)
Faecalibacterium_markers = FindMarkers(sample_fig, ident.1 = "Faecalibacterium", min.pct = 0.1)
Fusicatenibacter_markers = FindMarkers(sample_fig, ident.1 = "Fusicatenibacter", min.pct = 0.1)
CAG_markers = FindMarkers(sample_fig, ident.1 = "CAG-81", min.pct = 0.1)
Lachnospira_markers = FindMarkers(sample_fig, ident.1 = "Lachnospira", min.pct = 0.1)

Prevotella_markers_names = c("uvig-70434","uvig-574592","uvig-452546")
Phascolarctobacterium_markers_names = c("uvig-45508","uvig-435129","uvig-567159")
Dorea_markers_names = c("uvig-170449","uvig-243838","uvig-567223")
Clostridium_markers_names = c("ivig-2065","uvig-417742","uvig-311200")
Roseburia_markers_names = c("uvig-333967","uvig-440494","uvig-424669")
Faecalibacterium_markers_names = c("uvig-459064","uvig-50611","uvig-425896")
Fusicatenibacter_markers_names = c("uvig-159360","uvig-261210","uvig-494327")
CAG_markers_name = c("uvig-437397","uvig-14281","uvig-572239")
Lachnospira_markers_names = c("uvig-434731","uvig-151480","uvig-561074")

GPD_relation_prevotella <- c("uvig-429290","uvig-234321","uvig-139489")
GPD_relation_Dorea <- c("uvig-453418","uvig-148318","")
GPD_relation_Clostridium <- c("uvig-299095","ivig-1973","uvig-460329")
GPD_relation_Roseburia <- c("uvig-11401","uvig-55655","uvig-231605")
GPD_relation_Faecalibacterium <- c("uvig-253495","uvig-462454","uvig-229456")
GPD_relation_Fusicatenibacter <- c("ivig-1905","uvig-145092","uvig-146982")
GPD_relation_CAG <- c("uvig-242621","uvig-432243","uvig-245662")
GPD_relation_Lachnospira <- c("uvig-564459","uvig-324777","uvig-107148")

markers_name = c("uvig-70434","uvig-574592","uvig-452546","uvig-429290","uvig-234321","uvig-139489",
                "uvig-45508","uvig-435129","uvig-567159",
                "uvig-170449","uvig-243838","uvig-567223","uvig-453418","uvig-148318",
                "ivig-2065","uvig-417742","uvig-311200","uvig-299095","ivig-1973","uvig-460329",
                "uvig-333967","uvig-440494","uvig-424669","uvig-11401","uvig-55655","uvig-231605",
                "uvig-459064","uvig-50611","uvig-425896","uvig-253495","uvig-462454","uvig-229456",
                "uvig-159360","uvig-261210","uvig-494327","ivig-1905","uvig-145092","uvig-146982",
                "uvig-437397","uvig-14281","uvig-572239","uvig-242621","uvig-432243","uvig-245662",
                "uvig-434731","uvig-151480","uvig-561074","uvig-564459","uvig-324777","uvig-107148")
#color
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
pie(rep(6,n), col=sample(color, n))
col_vector
col_vector =c(wes_palette("Darjeeling1"), wes_palette("GrandBudapest1"), wes_palette("Cavalcanti1"), wes_palette("GrandBudapest2"), wes_palette("FantasticFox1"))
pal <- wes_palette("Zissou1", 10, type = "continuous")
pal2 <- wes_palette("Zissou1", 5, type = "continuous")
pal[3:10]

dotplot <- DotPlot(sample_fig, features = markers_name) +
           RotatedAxis() +
           #theme(axis.text.y = element_text(angle = 45, face="italic", hjust=1),
           theme(axis.text.y = element_text(face="italic", hjust=1,size=12),
           axis.text.x = element_text(size = 15)) + 
           scale_colour_gradientn(colours = pal) +
           theme(legend.position="right",
           legend.title = element_text(size=5),
           legend.text = element_text(size=5)) +
           #labs(title = "Genus related virus clusters", y = "", x = "") + 
           labs(title = "", y = "", x = "") +
           theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) +
           coord_flip() #x/y axis change

ggsave(filename = "dotplot_genus_related_virus_cluster.pdf",
      plot = dotplot,
      width = 6,
      height = 12,
      units = "in",
      dpi=500)


###############phage-related UMI count for each genus
Prevotella_subset = subset(sample_fig,idents = "Prevotella")
Phascolarctobacterium_A_subset = subset(sample_fig,idents = "Phascolarctobacterium_A")
Dorea_A_subset = subset(sample_fig,idents = "Dorea_A")
Clostridium_Q_subset = subset(sample_fig,idents = "Clostridium_Q")
Roseburia_subset = subset(sample_fig,idents = "Roseburia")
Faecalibacterium_subset = subset(sample_fig,idents = "Faecalibacterium")
Fusicatenibacter_subset = subset(sample_fig,idents = "Fusicatenibacter")
CAG_subset = subset(sample_fig,idents = "CAG-81")
Lachnospira_subset = subset(sample_fig,idents = "Lachnospira")

sum(Prevotella_subset$nCount_RNA)
sum(Phascolarctobacterium_A_subset$nCount_RNA)
sum(Dorea_A_subset$nCount_RNA)
sum(Clostridium_Q_subset$nCount_RNA)
sum(Roseburia_subset$nCount_RNA)
sum(Faecalibacterium_subset$nCount_RNA)
sum(Fusicatenibacter_subset$nCount_RNA)
sum(CAG_subset$nCount_RNA)
sum(Lachnospira_subset$nCount_RNA)

write.table(colnames(Prevotella_subset),"Prevotella_subset_bc.txt",row.names=F)
write.table(colnames(Phascolarctobacterium_A_subset),"Phascolarctobacterium_A_subset_bc.txt",row.names=F)
write.table(colnames(Dorea_A_subset),"Dorea_A_subset_bc.txt",row.names=F)
write.table(colnames(Clostridium_Q_subset),"Clostridium_Q_subset_bc.txt",row.names=F)
write.table(colnames(Roseburia_subset),"Roseburia_subset_bc.txt",row.names=F)
write.table(colnames(Faecalibacterium_subset),"Faecalibacterium_subset_bc.txt",row.names=F)
write.table(colnames(CAG_subset),"CAG_subset_bc.txt",row.names=F)
write.table(colnames(Fusicatenibacter_subset),"Fusicatenibacter_subset_bc.txt",row.names=F)
write.table(colnames(Lachnospira_subset),"Lachnospira_subset_bc.txt",row.names=F)


###phage-related component rate for each barcode
phage_rate<- read.csv("/public2/labmember/wyc_rawdata/A197/supplement_data/phage_test/nonrrna_star_test/result/A197_bc_phage_rate.txt",header = F,sep = "\t") 
colnames(phage_rate) = c("BC","map","sum","rate")
BC <- colnames(sample_fig)
BC<-as.data.frame(BC)
phage_rate_info <- merge(BC, phage_rate,all.x=TRUE)
phage_rate_info[is.na(phage_rate_info)] = 0
sample_fig@meta.data$phage_rate  <- phage_rate_info[,4]

plot <- VlnPlot(sample_fig, features = "phage_rate", col=cell_type_cols,pt.size = 0)+ 
        NoLegend() +
        geom_boxplot(width=.1, col="black", lwd = 1.2, fatten = 1.5, outlier.shape = NA) + 
        theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 30),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=30)) + 
        theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

ggsave(filename = "phage_rate.pdf",
      plot = plot,
      width = 15,
      height = 7,
      units = "in",
      dpi=500)

