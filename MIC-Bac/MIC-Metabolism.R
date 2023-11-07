library(KEGGREST)
library(ggplot2)
library(GSVA)
library(GSEABase)
library(pheatmap)
library(ComplexHeatmap)
library(stringr)
## merge select genes and uhgg_genome
genome_uhgg <- read.table("uhgg_genome.txt", sep = "\t",header = T,stringsAsFactors = FALSE,quote = "")
colnames(genome_uhgg) <- c("gene","function")
gene_num <- select_gene_list
merge(gene_num,genome_uhgg, by = "gene") -> list_gene

## get the product, gene name and KEGGnumber of select genes 
get_function <- list_gene
get_function_l <- strsplit(get_function$function.,split = ";")

gene_f <- matrix(data = NA,nrow = length(get_function_l),ncol = 3)
for (i in 1:length(get_function_l)){
  for (si in 1:length(get_function_l[[i]])){
    if(str_detect(get_function_l[[i]][si],"KEGG")){gene_f[i,1] <- get_function_l[[i]][si]}
    if(str_detect(get_function_l[[i]][si],"Name")){gene_f[i,2] <- get_function_l[[i]][si]}
    if(str_detect(get_function_l[[i]][si],"product")){gene_f[i,3] <- get_function_l[[i]][si]}
  }
}
gene_f <- as.data.frame(gene_f)
## get KEGG number, gene name and product 
colnames(gene_f) <- c("KEGGnumber","genename",'product')
get_function <- cbind(get_function,gene_f)
get_function <- na.omit(get_function)
get_function <- get_function[!duplicated(get_function$gene),]
## discard ribosomal genes
grep("ribosomal",get_function$product) -> id
get_function <- get_function[-c(id),]
## discard genes without KEGG number
get_function <- get_function[-c(which(get_function$KEGGnumber == 'KEGG=-')),]
get_function_m <- strsplit(get_function$KEGGnumber,split = ",")
get_single_number <- matrix(data = NA,nrow = length(get_function_m),ncol = 1)
for (i in 1:length(get_function_m)){
  get_single_number[i,1] <- get_function_m[[i]][1]
}
get_single_number <- as.data.frame(get_single_number)
colnames(get_single_number) <- "KEGGnumber"
get_function$KEGGnumber <- get_single_number$KEGGnumber
get_function <- get_function[!duplicated(get_function$KEGGnumber),]
get_function$KEGGnumber <- gsub("KEGG=",'',get_function$KEGGnumber)
get_function$genename <- gsub("Name=",'',get_function$genename)
get_function$product <- gsub('product=','',get_function$product)
##get pathway
rownames(get_function) <- get_function$gene

kegg_path <- matrix(data = NA,nrow = nrow(get_function),ncol = 20)
for (i in 1:nrow(get_function)){
  if(i%%10==0){print(i)}
  query <- keggGet(c(get_function$KEGGnumber)[i])
  for (si in 1:length(query)){
    names(query[[1]]$PATHWAY) -> qi
    if(is.null(qi)){kegg_path[i,1] <- "nopath"}else{kegg_path[i,c(1:length(query[[1]]$PATHWAY))] <- qi }
  }
}
kegg_path_n <- as.data.frame(kegg_path)
kegg_path_n <- cbind(get_function,kegg_path_n)
##discard no path
kegg_path_n <- kegg_path_n[-c(which(kegg_path_n$V1=="nopath")),]

cc <- c()
for (i in 1:nrow(kegg_path_selg)){
  for (s in 1:(which(is.na(kegg_path_selg[i,]))[1]-1)){
    kegg_path_selg[i,s] -> ci
    cc <- rbind(cc,ci)
  }
}
cc <- unique(cc)

list_path <- matrix(data = NA,nrow = length(cc),ncol = nrow(kegg_path_selg))
rownames(list_path) <- cc
for (i in 1:nrow(list_path)){
  for (si in 1:nrow(kegg_path_selg)){
    for (q in 1:(which(is.na(kegg_path_selg[si,]))[1]-1)){
      if(kegg_path_selg[si,q]==rownames(list_path)[i]){list_path[i,si] <- rownames(kegg_path_selg)[si]}
    }
  }
}
list_path <- as.data.frame(list_path)

## get KEGG pathway annotation
kegg_pathway_anno <- c()
for (i in 1:nrow(list_path)){
  if(i%%10==0){print(i)}
  path1 <- keggGet(c(rownames(list_path)[i])) 
  path1[[1]]$PATHWAY_MAP -> path2
  kegg_pathway_anno <- rbind(kegg_pathway_anno,path2)
}
list_path$pathwayanno <- kegg_pathway_anno[,1]

gmtFile <- "metabolism.gmt"
countexp2 <- as.data.frame(scale.data)
rownames(countexp2) <- gsub("-","_",rownames(countexp2))
geneSets <- getGmt(gmtFile) #signature read
gsva_es <- gsva(as.matrix(countexp2), geneSets, method=c("ssgsea"), kcdf=c("Poisson")) 
signature_exp<-data.frame(gsva_es)

selg_path_t <- as.data.frame(t(signature_exp))
avg_df = aggregate(selg_path_t[,1:ncol(selg_path_t)],list(metadata$orig.ident),mean)
rownames(avg_df) <- avg_df$Group.1
avg_df <- avg_df[,-1]
avg_df <- t(avg_df)
## color 
heatmap_color=RColorBrewer::brewer.pal(name = "RdBu",n=11)
pal=rev(colorRampPalette(heatmap_color)(100))
## draw
p <- pheatmap(t(avg_df), 
              color = pal,
              border_color = NA,scale = "row",
              #labels_col="",
              cluster_cols = F,show_colnames = T,cluster_rows = F,treeheight_row = 0,
              fontsize=10, fontsize_row=12,fontsize_col=10) #set color

p

