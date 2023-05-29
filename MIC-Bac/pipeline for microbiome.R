library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
library(tidyverse)


##filt count matrix 
data_A197 <- read.csv(file = "/public3/project_users/shenyf/brain/qu/singlecell/genus_info/A197/A197_multi_cell_matrix.txt", sep = "\t", header = T, row.names = 1)
genus_info_A197 <- read.csv(file = "/public2/labmember/wyc_rawdata/A197/supplement_data/phage_test/star_test/bc_info.txt" , header = F, sep = "\t")

cluster_pipeline <- function(input_data=NA, bacteria_info = NA, min_cells=30, min_features=15,
                             sp_rate_threshold = 0.03, pc.num=1:30, resolution_num = 0.5, if_nr = FALSE) {
  ##create Seurat object
  sample <- CreateSeuratObject(counts = input_data, project = "A197",min.cells = min_cells, min.features = min_features)
  ##get annotation
  BC <- colnames(sample)
  BC<-as.data.frame(BC)
  colnames(bacteria_info) <- c("BC","index")
  BC_info <- merge(BC, bacteria_info,all.x=TRUE)
  sample@meta.data$genus_info <- BC_info[,2]
  genus_num_data = table(sample@meta.data$genus_info)
  
  ##filt low abundance bacteria
  num_threshold = dim(sample)[2] * sp_rate_threshold
  for(i in 1:length(genus_num_data)){
    if(genus_num_data[[i]]<num_threshold){
      index = names(genus_num_data[i])
      sample <- subset(sample,subset = genus_info != index)
    }
  }
  ##reduction
  sample <- NormalizeData(sample, normalization.method = "LogNormalize", scale.factor = 10000)
  sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000)
  sample <- ScaleData(sample,features = VariableFeatures(sample))
  sample <- RunPCA(sample,features = VariableFeatures(sample))
  sample <- FindNeighbors(sample,dims = pc.num)
  sample <- FindClusters(sample, resolution = resolution_num)
  sample <- RunUMAP(sample, dims = pc.num)
  return(sample)
}
sample_newinfo <- cluster_pipeline(input_data = data_A197, bacteria_info = genus_info_A197, min_cells = 50, min_features = 30) 

saveRDS(sample_newinfo,"./sample.RDS")

####inload .RDS data
sample <- readRDS("C:\\Users\\86153\\Desktop\\qu\\single_micro\\barcode_genus\\A197\\qqh\\sample.RDS")

#color set
cell_type_cols <- c(brewer.pal(9, "Set1"), "#FF34B3","#BC8F8F","#20B2AA","#00F5FF",
                    "#FFA500","#ADFF2F","#FF6A6A","#7FFFD4", "#AB82FF","#90EE90",
                    "#00CD00","#008B8B","#6495ED","#FFC1C1","#CD5C5C","#8B008B",
                    "#FF3030", "#7CFC00","#000000","#708090")

sample_meta <- as.data.frame(sample@meta.data)
##change name
sample_meta$genus_info[which(sample_meta$genus_info=="Clostridium_Q")] <- "Clostridium"
sample_meta$genus_info[which(sample_meta$genus_info=="Dorea_A")] <- "Dorea"
sample_meta$genus_info[which(sample_meta$genus_info=="Phascolarctobacterium_A")] <- "Phascolarctobacterium"
##add new name in metadata
sample@meta.data$genus_info_new <- sample_meta$genus_info
Idents(sample) <- sample$genus_info_new

##get species annotation
meta_species <- read.table("C:\\Users\\86153\\Desktop\\qu\\single_micro\\barcode_species\\A197\\qqh\\A197.barcode_count.txt",sep = "\t",row.names = 1)
meta_sample <- as.data.frame(sample@meta.data)
ip <- intersect(rownames(meta_sample),rownames(meta_species))
meta_species <- meta_species[c(ip),]

meta_species$V6->spn
strsplit(spn,split="\\|")->xps
kwx<-c()
for(i in 1:length(spn)){
  
  kw<-xps[[i]][length(xps[[i]])-1]
  kwx<-c(kwx,kw)
  
}

ctxn<-data.frame(kwx,meta_species)

meta_species$V6->spn
strsplit(spn,split="\\|")->xps
kwx<-c()
for(i in 1:length(spn)){
  
  kw<-xps[[i]][length(xps[[i]])]
  
  kwx<-c(kwx,kw)
  
}
ctxx<-data.frame(kwx,ctxn)
colnames(ctxx)[1] <- "species"
colnames(ctxx)[2] <- "genus"
ctxx <- as.data.frame(ctxx[,c(1,2)])
ctxx$species <- gsub("s_","",ctxx$species)
ctxx$genus <- gsub("g_","",ctxx$genus)

#meta Phascolarctobacterium_A
meta_pha <- ctxx
meta_pha$species[which(ctxx$genus!="Phascolarctobacterium_A")] <- "others"

#meta prevotella
meta_pre <- ctxx
meta_pre$species[which(ctxx$genus!="Prevotella")] <- "others"
meta_pre_new <- meta_pre
which(meta_pre_new$species %in% c("Prevotella hominis","Prevotella sp000434975","Prevotella sp000436035","Prevotella sp900313215","Prevotella sp900544825",
                                  "Prevotella sp900546535","Prevotella sp900555035","Prevotella sp900556795","Prevotella sp900557255","Prevotella sp900765465","Prevotella stercorea")) -> id

meta_pre_new$species[id] <- "Prevotella_other"

#meta Roseburia
meta_ros <- ctxx
meta_ros$species[which(ctxx$genus!="Roseburia")] <- "others"

#add metadata
sample@meta.data$species_info_pha <- meta_pha[,1]
sample@meta.data$species_info_ros <- meta_ros[,1]
sample@meta.data$species_info_pre_new <- meta_pre_new[,1]


#####Phascolarctobacterium_A succinatutens analyse
sample_pha <- sample[,sample@meta.data$species_info_pha %in% c("Phascolarctobacterium_A succinatutens")]
#reduction
sample_pha <- NormalizeData(sample_pha, normalization.method = "LogNormalize", scale.factor = 10000)
sample_pha <- FindVariableFeatures(sample_pha, selection.method = "vst", nfeatures = 2000)
sample_pha <- ScaleData(sample_pha,features = VariableFeatures(sample_pha))
sample_pha <- RunPCA(sample_pha,features = VariableFeatures(sample_pha))
sample_pha <- FindNeighbors(sample_pha,dims = 1:30)
sample_pha <- FindClusters(sample_pha, resolution = 0.5)
sample_pha <- RunUMAP(sample_pha, dims = 1:30)

cluster2celltype <- c("0"="0", "1"="1","2"="2","3"="0")
sample_pha[['anno1_CellType']] = unname(cluster2celltype[sample_pha@meta.data$seurat_clusters])
Idents(sample_pha) <- sample_pha$anno1_CellType
levels(sample_pha) <- c("0","1","2")

#clusters in species  Phascolarctobacterium_A succinatutens
DimPlot(sample_pha, reduction = "umap", group.by = 'anno1_CellType', cols = cell_type_cols) +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(axis.text.y = element_blank(),
        text = element_text(size = 20),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 20),
        axis.ticks.x = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

#get marker gene
sample_pha_markers <- FindAllMarkers(sample_pha, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1)
#get gene expression UMAP
FeaturePlot(sample_pha, features = c("MGYG000001365-01624"))
##get gene expression dot plot 
markers <- c("MGYG000001365-00002","MGYG000001365-00617","MGYG000001365-01624","MGYG000001365-01277","MGYG000001365-01301","MGYG000001365-00508","MGYG000001365-01761","MGYG000001365-01923","MGYG000001365-00149",
             "MGYG000000591-01347","MGYG000001365-00414","MGYG000001365-00415","MGYG000001365-01214")
sample_meta_t <- sample_pha@meta.data
sample_meta_t$CB <- rownames(sample_meta_t)
sample_pha$CB <- sample_meta_t$CB
bubble.df=as.matrix(sample_pha[["RNA"]]@data[markers,])
bubble.df=t(bubble.df)
bubble.df=as.data.frame(scale(bubble.df))
colnames(bubble.df) <- c("ISClte1","IS663 ","ISL7","mdtK","mdtB","mdtG","mdtC","rcsC_2","rcsC_1",
                         "scpA_3","mutB","scpA_1","pccB")
bubble.df$CB=rownames(bubble.df)
bubble.df=merge(bubble.df,sample_pha@meta.data[,c("CB","anno1_CellType")],by = "CB")
bubble.df$CB=NULL
markers <- c("ISClte1","IS663 ","ISL7","mdtK","mdtB","mdtG","mdtC","rcsC_2","rcsC_1",
             "scpA_3","mutB","scpA_1","pccB")

celltype_v=c()
gene_v=c()
mean_v=c()
ratio_v=c()
for (i in unique(bubble.df$anno1_CellType)) {
  bubble.df_small=bubble.df%>%filter(anno1_CellType==i)
  for (j in markers) {
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
plotdf$gene=factor(plotdf$gene,levels = rev(as.character(markers)))
plotdf$exp=ifelse(plotdf$exp>3,3,plotdf$exp)
plotdf%>%ggplot(aes(x=celltype,y=gene,size=ratio,color=exp))+geom_point()+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11])))+
  scale_size_continuous(limits = c(0,1))+theme_bw()+scale_fill_continuous(breaks=c(0,0.5,1))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )

##correlation analyse
###IS family matrix
markers <- c("MGYG000001365-00002","MGYG000001365-00617","MGYG000001365-01624")
sample_pha_is=as.matrix(sample_pha[["RNA"]]@data[markers,])
sample_pha_is=t(sample_pha_is)
sample_pha_is=as.data.frame(scale(sample_pha_is))
colnames(sample_pha_is) <- c("IS3","IS1182","IS30")
##mdt matrix
markers <- c("MGYG000001365-01277","MGYG000001365-01301","MGYG000001365-00508","MGYG000001365-01761")
sample_pha_mdt <- as.matrix(sample_pha[["RNA"]]@data[markers,])
sample_pha_mdt=t(sample_pha_mdt)
sample_pha_mdt=as.data.frame(scale(sample_pha_mdt))
colnames(sample_pha_mdt) <- c("mdtK","mdtB","mdtG","mdtC")

species_list<-c()
meta_list<-c()
p_list<-c()
r_list<-c()
comat_p<- matrix(nrow=ncol(sample_pha_mdt),ncol=(ncol(sample_pha_is)))
comat_r<- matrix(nrow=ncol(sample_pha_mdt),ncol=(ncol(sample_pha_is)))
for(si in 1:ncol(sample_pha_mdt)){
  if(si %% 10 == 0 ){print(si)}
  si_name<-colnames(sample_pha_mdt)[si]
  si_data<-sample_pha_mdt[,si]
  for(mi in 1:ncol(sample_pha_is)){
    mi_name<-colnames(sample_pha_is)[mi]
    mi_data<-sample_pha_is[,mi]
    cor_result<-cor.test(as.numeric(si_data), as.numeric(mi_data),  method = "pearson")
    p_result<-cor_result$p.value
    r_result<-cor_result$estimate
    comat_p[si,mi]<-p_result
    comat_r[si,mi]<-r_result
    if(p_result<=0.05){
      species_list<-c(species_list,si_name)
      meta_list<-c(meta_list,mi_name)
      p_list<-c(p_list,p_result)
      r_list<-c(r_list,r_result)
    }
  }
}

sig_data<-data.frame(species=species_list,meta=meta_list,pvalue=p_list,cor=r_list)

rownames(comat_p)<-colnames(sample_pha_mdt)
colnames(comat_p)<-colnames(sample_pha_is)[1:ncol(sample_pha_is)]
rownames(comat_r)<-colnames(sample_pha_mdt)
colnames(comat_r)<-colnames(sample_pha_is)[1:ncol(sample_pha_is)]
comat_r_min <- comat_r
comat_p_min <- comat_p
colnames(comat_p_min) <- c("ISClte1","IS663 ","ISL7")
colnames(comat_r_min) <- c("ISClte1","IS663 ","ISL7")
library(corrplot)
COLOR01<- colorRampPalette(c("#5081ff","#638dff","#b1c6ff","#c3d5ff" ,"#ffffff", "#f7cccc","#f1aeac","#eb8b8a","#e66a68"))(50)
corrplot(comat_r_min,is.corr = T,
         type="full",
         p.mat=comat_p_min,
         insig = "label_sig",
         sig.level = c(0.001,0.01,0.05),
         pch.cex = 1.6,pch.col="black",
         #diag = T,
         tl.srt =50,
         tl.col = "black",
         #family="serif",
         col = COLOR01,
         tl.cex = 1,title = "",cl.lim=c(-1,1),
         cl.pos="b",cl.length=3,cl.ratio=0.1,
         mar=c(0, 7, 1, 7),cl.cex = 0.8, method = "color"
)

###get genus UMAP
sample_new <- sample[,sample@meta.data$RNA_snn_res.0.5 %in% c(0,1,4,5,6,7,8,9,10,11,12)]

cluster2celltype <- c("0"="0", "1"="1","2"="0","3"="0","4"="2", "5"="3", "6"="4", "7"="5", "8"="6",
                      "9"="7","10"="8", "11"="9", "12"="10", "13"="11",
                      "14"="12")

sample_new[['anno1_CellType']] = unname(cluster2celltype[sample_new@meta.data$seurat_clusters])
sample_new$genus_info <- factor(levels = c("Prevotella","Phascolarctobacterium_A","Clostridium_Q","Dorea_A","Roseburia","Lachnospira","Faecalibacterium","CAG-81","Fusicatenibacter"),sample_new$genus_info)
DimPlot(sample_new, reduction = "umap",cols= cell_type_cols,pt.size = 0.005 ,group.by = 'genus_info') +
  #NoLegend() +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(axis.text.y = element_blank(),
        text = element_text(size = 15),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 15),
        axis.ticks.x = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))

##genus Prevotella analyse
#get marker gene
sample_pre <- sample_new[,sample_new@meta.data$species_info_pre_new %in% c("Prevotella copri_A","Prevotella sp900767615","Prevotella_other")]
Idents(sample_pre) <- sample_pre$species_info_pre_new
sample_pre_markers <- FindAllMarkers(sample_pre, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1)
#get expression dot plot
markers <- c("MGYG000002603-00091","MGYG000002603-01642","MGYG000002603-01855","MGYG000002960-00547","MGYG000002960-00933","MGYG000003347-00303")
sample_meta_t <- sample_pre@meta.data
sample_meta_t$CB <- rownames(sample_meta_t)
sample_pre$CB <- sample_meta_t$CB
bubble.df=as.matrix(sample_pre[["RNA"]]@data[markers,])
bubble.df=t(bubble.df)
bubble.df=as.data.frame(scale(bubble.df))
colnames(bubble.df) <- c("rcsC_1","rcsC_5","rcsC_7","rcsC_2","rcsC_5","rcsC_3")
bubble.df <- bubble.df[,-2]
bubble.df$CB=rownames(bubble.df)
bubble.df=merge(bubble.df,sample_pre@meta.data[,c("CB","species_info_pre_new")],by = "CB")
markers <- c("rcsC_1","rcsC_7","rcsC_2","rcsC_5","rcsC_3")
celltype_v=c()
gene_v=c()
mean_v=c()
ratio_v=c()
for (i in unique(bubble.df$species_info_pre_new)) {
  bubble.df_small=bubble.df%>%filter(species_info_pre_new==i)
  for (j in markers) {
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
plotdf$gene=factor(plotdf$gene,levels = rev(as.character(markers)))
plotdf$exp=ifelse(plotdf$exp>3,3,plotdf$exp)
plotdf%>%ggplot(aes(x=celltype,y=gene,size=ratio,color=exp))+geom_point()+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11])))+
  scale_size_continuous(limits = c(0,1))+theme_bw()+scale_fill_continuous(breaks=c(0,0.5,1))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45),axis.text.x = element_text(face = "italic")
  )

##genus Roseburia analyse
sample_ros <- sample_new[,sample_new@meta.data$genus_info %in% c("Roseburia")]
Idents(sample_ros) <- sample_ros$species_info_ros
levels(sample_ros)
#get marker gene
sample_ros_markers <- FindAllMarkers(sample_ros, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1)
sample_meta_t <- sample_ros@meta.data
sample_meta_t$CB <- rownames(sample_meta_t)
sample_ros$CB <- sample_meta_t$CB

markers <- c("MGYG000002517-00135","MGYG000000076-00769","MGYG000000076-01417")
bubble.df=as.matrix(sample_ros[["RNA"]]@data[markers,])
bubble.df=t(bubble.df)
bubble.df=as.data.frame(scale(bubble.df))
colnames(bubble.df) <- c("hag_1","sasA_3","mdtC")
bubble.df$CB=rownames(bubble.df)
bubble.df=merge(bubble.df,sample_ros@meta.data[,c("CB","species_info_ros")],by = "CB")
bubble.df$CB=NULL
markers <- c("hag_1","sasA_3","mdtC")

celltype_v=c()
gene_v=c()
mean_v=c()
ratio_v=c()
for (i in unique(bubble.df$species_info_ros)) {
  bubble.df_small=bubble.df%>%filter(species_info_ros==i)
  for (j in markers) {
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
plotdf$gene=factor(plotdf$gene,levels = rev(as.character(markers)))
plotdf$exp=ifelse(plotdf$exp>3,3,plotdf$exp)
plotdf%>%ggplot(aes(x=celltype,y=gene,size=ratio,color=exp))+geom_point()+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11])))+
  scale_size_continuous(limits = c(0,1))+theme_bw()+scale_fill_continuous(breaks=c(0,0.5,1))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45),
    axis.text.x = element_text(face = "italic")
  )
#UMAP colored by species
sample_draw <- sample_new
sample_meta <- sample_draw@meta.data
which(sample_meta$genus_info=="Prevotella") -> id
sample_meta$species_info_ros[id] <- sample_meta$species_info_pre_new[id]
sample_draw$species_combine <- sample_meta$species_info_ros
sample_draw$species_combine <- factor(level = c("Prevotella copri_A","Prevotella sp900767615","Prevotella_other","Roseburia hominis","Roseburia intestinalis","Roseburia sp900552665","others"),sample_draw$species_combine)
DimPlot(sample_draw, reduction = "umap",cols= c("#377EB8","#4DAF4A","#984EA3","#E41A1C","#FF7F00","#FFFF33","#999999"),group.by = 'species_combine') +
  #NoLegend() +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(axis.text.y = element_blank(),
        text = element_text(size = 15),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 15),
        axis.ticks.x = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))








