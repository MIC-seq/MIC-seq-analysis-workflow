library(vegan)
library(reshape2)
library(paletteer)
library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(GSVA)
library(GSEABase)
library(gplots)
library(VennDiagram)
library(monocle)
## color
cell_type_cols <- c(brewer.pal(9, "Set1"), "#FF34B3","#BC8F8F","#20B2AA","#00F5FF",
                    "#FFA500","#ADFF2F","#FF6A6A","#7FFFD4", "#AB82FF","#90EE90",
                    "#00CD00","#008B8B","#6495ED","#FFC1C1","#CD5C5C","#8B008B",
                    "#FF3030", "#7CFC00","#000000","#708090")
heatmap_color=RColorBrewer::brewer.pal(name = "RdBu",n=11)
pal=rev(colorRampPalette(heatmap_color)(100))

## load a demo of three donor dataset
combine_sample <- readRDS('combine_dataset_demo.rds')

## clustering and re-clustering ##########
pc.num <- 1:30
clustering_analysis <- function(sample=NA,resolution_num=NA){
  sample <- NormalizeData(sample, normalization.method = "LogNormalize", scale.factor = 10000)
  sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 3000)
  sample <- ScaleData(sample,features = VariableFeatures(sample))
  sample <- RunPCA(sample,features = VariableFeatures(sample))
  sample <- FindNeighbors(sample,dims = pc.num, reduction = "pca")
  sample <- FindClusters(sample, resolution = resolution_num)
  sample <- RunUMAP(sample, dims = pc.num, reduction = "pca")
  return(sample)
}

combine_sample <- clustering_analysis(combine_sample,resolution_num = 0.5)
DimPlot(combine_sample, reduction = "umap",label.size = 6,label = F,cols = cell_type_cols,
        group.by = 'species_info') +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(axis.text.y = element_blank(),
        text = element_text(size = 20),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

## using for draw marker gene dotplot #########
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

####### MIC-Metabolism ############
#gsva function 
#1.load gmt file
#use M. funiformis as an example

gmtFile_mega <- "gmt_mega.gmt"
geneSets_mega_total <- getGmt(gmtFile_mega) #signature read
mega_count_matrix <- readRDS('mega_count_matrix.rds')


#1.Raw MES generation # MIC-Metabolism  supports different methods: ssGSEA, AUCell , VISION and GSVA 

raw_mes <- function(geneset_input = NA,count_matrix = NA,method = "ssGSEA"){
       if (method == "ssGSEA") {
           library(GSVA)
           library(GSEABase)
           gsva_es <- gsva(as.matrix(count_matrix), geneset_input,method=c("ssgsea"), kcdf=c("Poisson")) 
           gsva_es <- as.data.frame(gsva_es)
         }
       if (method == "GSVA") {
           library(GSVA)
           library(GSEABase)
           gsva_es <- gsva(as.matrix(count_matrix), geneset_input, method=c("gsva"), kcdf=c("Poisson"))
           gsva_es <- as.data.frame(gsva_es)
         }
       if (method == "AUCell") {
           library(AUCell)
           library(GSEABase)
           cells_rankings <- AUCell_buildRankings(as.matrix(count_matrix), plotStats=F)
           cells_AUC <- AUCell_calcAUC(geneset_input, cells_rankings)
           gsva_es <- data.frame(getAUC(cells_AUC))
         }
       if (method == "VISION") {
           library(VISION)
           n.umi <- colSums(count_matrix)
           scaled_counts <- t(t(count_matrix) / n.umi) * median(n.umi)
           vis <- Vision(scaled_counts, signatures = geneset_input)
           vis <- analyze(vis)
           gsva_es<-data.frame(t(vis@SigScores))
         }
       return(gsva_es)
     }
gsva_es <- raw_mes(geneset_input = geneSets_mega_total,count_matrix = mega_count_matrix,method = "ssGSEA")

#2.perform permutation
get_gene_count <- function(geneset_input = NA){
  gene_count <- matrix(data = NA,nrow = length(geneset_input),ncol = 2)
  for (i in 1:length(geneset_input)){
    length(c(geneset_input[[i]]@geneIds)) -> gene_count[i,1]
    gene_count[i,2] <- geneset_input[[i]]@setName
  }
  
  gene_count <- as.data.frame(gene_count)
  colnames(gene_count) <- c('gene_number','pathway_name')
  return(gene_count)
}
gmt_gene_count <- get_gene_count(geneset_input = geneSets_mega_total)

get_random_geneset <- function(pathway_gsva = NA,gene_count = NA,count_matrix = NA){
  total_pathway_random <- c()
  for (q in 1:nrow(pathway_gsva)){
    if (q %%10==0){print(q)}
    gene_select_number <- gene_count$gene_number[which(gene_count$pathway_name==rownames(pathway_gsva)[q])]
    gene_select_number <- as.numeric(gene_select_number)
    pathway_random <- c()
    for (i in 1:1000){
      sample(rownames(count_matrix),gene_select_number) -> gene_select
      pathway_random[[i]] <- c(gene_select)
      names(pathway_random)[[i]] <- paste('random',i,sep = "_")
    }
    pathway_random -> total_pathway_random[[rownames(pathway_gsva)[q]]]
  }
  return(total_pathway_random)
}

random_set <- get_random_geneset(pathway_gsva = gsva_es,gene_count = gmt_gene_count,count_matrix = mega_count_matrix)

get_permutation <- function(pathway_gsva = NA,count_matrix = NA,total_pathway_random = NA){
  first_permutation <- c()
  for (q in 1:nrow(pathway_gsva)){
    if (q %%10==0){print(q)}
    gene_gsva <- gsva(as.matrix(count_matrix), total_pathway_random[[rownames(pathway_gsva)[q]]],method=c("ssgsea"), kcdf=c("Poisson"))
    first_permutation[[rownames(pathway_gsva)[q]]] <- gene_gsva
  }
  return(first_permutation)
}
first_permutation <- get_permutation(pathway_gsva = gsva_es,count_matrix = mega_count_matrix,total_pathway_random = random_set)

#3.scale
# for each species
gsva_es_mean <- as.data.frame(apply(gsva_es,1,mean))
permutation_first_result <- matrix(data = NA,nrow = length(first_permutation),ncol = 1)
get_zcore <- matrix(data = NA,nrow = length(first_permutation),ncol = 1)
for (i in 1:length(first_permutation)){
  rank_data <- as.data.frame(first_permutation[[i]])
  apply(rank_data,1,mean) -> rank_data
  rank_data <- as.data.frame(rank_data)
  rank_data <- rbind(gsva_es_mean[i,1],rank_data)
  get_zcore_n <- scale(rank_data$rank_data)
  get_zcore[i,1] <- get_zcore_n[1,1]
  rank_test <- rank(-(rank_data[c(1:1001),]))
  for (q in 1:length(rank_test)){
    if(rank_test[q]>1){rank_test[q] <- rank_test[q]-1}
  }
  rank_test <- rank_test/1000
  rank_test[1] -> permutation_first_result[i,1]
}
permutation_first_result <- as.data.frame(permutation_first_result)
rownames(permutation_first_result) <- rownames(gsva_es_mean)
get_zcore <- as.data.frame(get_zcore)
rownames(get_zcore) <- rownames(gsva_es_mean)


######### pseudotime analysis #########
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
diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~sct_cluster_new", cores = 1) 
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
















