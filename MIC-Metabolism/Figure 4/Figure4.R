#### Figure 4

### MIC-matabolism

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

## Fig.4E
p <- ggplot(gsva_score,aes(x=UMAP_1,y=UMAP_2))+
  
  geom_point(data = mydata,aes(x=UMAP_1,y=UMAP_2,
                               
                               color=pyruvate_score),size=0.1)+
  
  scale_colour_gradientn('selg',limits=c(0,0.2),colors = c('#EAEBEB','yellow','red','darkred'))+
  
  theme_bw()+theme(panel.grid.major=element_line(colour=NA),
                   
                   panel.background = element_rect(fill = "transparent",colour = NA),
                   
                   plot.background = element_rect(fill = "transparent",colour = NA),
                   
                   panel.grid.minor = element_blank())
p

ggplot(pyruvate_gsva,aes(x=species, y=Pyruvate_metabolism,color=species,fill=species,alpha=0.6))+ #scale_y_log10(limits = c(1, 100))+ 
  geom_boxplot(color = c(cell_type_cols[1:5]),size=1,width = 0.5) + theme_bw() + 
  
  
  theme(panel.grid = element_blank(),
        
        panel.background = element_blank(),
        
        axis.line = element_line(),
        
        axis.text.x = element_text(size=12,color="black", face= "bold",angle=60,hjust = 1),
        
        
        axis.text.y = element_text(size=12, color="black", face= "bold"),
        
        plot.title = element_text(size=12,color="black", face= "bold",hjust=0))+
  
  #geom_jitter(width = 0.2,show.legend =F,size=1.2) +
  
  scale_fill_manual(values=c(cell_type_cols[1:5]))+
  
  labs(x="",y="",title="")+#ylim(0.1,0.3)+
  
  theme(strip.text.x = element_text(size = 12, colour = "black"))+
  
  theme(strip.background.x = element_rect(fill = "lightgrey")) 

### Fig. 4F
p <- ggplot(gsva_score,aes(x=UMAP_1,y=UMAP_2))+
  
  geom_point(data = mydata,aes(x=UMAP_1,y=UMAP_2,
                               
                               color=pro_score),size=0.1)+
  
  scale_colour_gradientn('selg',limits=c(0,0.2),colors = c('#EAEBEB','yellow','red','darkred'))+
  
  theme_bw()+theme(panel.grid.major=element_line(colour=NA),
                   
                   panel.background = element_rect(fill = "transparent",colour = NA),
                   
                   plot.background = element_rect(fill = "transparent",colour = NA),
                   
                   panel.grid.minor = element_blank())
p

ggplot(Propanoate_gsva,aes(x=species, y=Propanoate_metabolism,color=species,fill=species,alpha=0.6))+ #scale_y_log10(limits = c(1, 100))+ 
  geom_boxplot(color = c(cell_type_cols[1:5]),size=1,width = 0.5) + theme_bw() + 
  
  
  theme(panel.grid = element_blank(),
        
        panel.background = element_blank(),
        
        axis.line = element_line(),
        
        axis.text.x = element_text(size=12,color="black", face= "bold",angle=60,hjust = 1),
        
        
        axis.text.y = element_text(size=12, color="black", face= "bold"),
        
        plot.title = element_text(size=12,color="black", face= "bold",hjust=0))+
  
  #geom_jitter(width = 0.2,show.legend =F,size=1.2) +
  
  scale_fill_manual(values=c(cell_type_cols[1:5]))+
  
  labs(x="",y="",title="")+#ylim(0.1,0.3)+
  
  theme(strip.text.x = element_text(size = 12, colour = "black"))+
  
  theme(strip.background.x = element_rect(fill = "lightgrey")) 


### Fig. 4G
selg_pathway <- c('Citrate cycle TCA cycle',"Pentose phosphate pathway","Glycolysis",'Pyruvate metabolism')

b_zscore_Et_b_selg <- b_zscore_Et_b[c(selg_pathway),]
b_zscore_Et_F_selg <- b_zscore_Et_F[c(selg_pathway),]


b_zscore_Et_b_selg <- as.data.frame(t(b_zscore_Et_b_selg))

b_zscore_Et_F_selg <- as.data.frame(t(b_zscore_Et_F_selg))


b_zscore_Et_b_selg$donor <- "ET-B"
b_zscore_Et_F_selg$donor <- "ET-F"

B_zscore_selg <- rbind(b_zscore_Et_b_selg,b_zscore_Et_F_selg)

B_zscore_selg <- melt(B_zscore_selg)

B_zscore_selg$variable <- gsub(' ','_',B_zscore_selg$variable)

B_zscore_selg$variable <- gsub('/','',B_zscore_selg$variable)

B_zscore_selg$variable <- gsub('[(]','',B_zscore_selg$variable)

B_zscore_selg$variable <- gsub('[)]','',B_zscore_selg$variable)

ggplot(B_zscore_selg, aes(x = variable, y = value))+ 
  
  geom_boxplot(aes(fill = donor),position=position_dodge(0.6),width=0.6,outlier.color = NA)+
  
  scale_fill_manual(values = c(cell_type_cols[1:7]))+theme_bw()+
  
  theme(panel.grid = element_blank(),
        
        panel.background = element_blank(),
        
        axis.line = element_line(),
        axis.text.x = element_text(size=12,color="black", face= "bold",angle=60,hjust = 1),
        
        axis.text.y = element_text(size=12, color="black", face= "bold"),
        
        plot.title = element_text(size=12,color="black", face= "bold",hjust=0))+
  
  stat_compare_means(aes(group = donor),
                     
                     label="p.signif",
                     
                     show.legend = F)


### Fig. 4H
selg_pathway <- c('Arginine biosynthesis',"Protein export","Bacterial secretion system",'Methane metabolism')

mega_zscore_d1_selg <- mega_zscore_d1[c(selg_pathway),]

mega_zscore_d3_selg <- mega_zscore_d3[c(selg_pathway),]

mega_zscore_d1_selg <- as.data.frame(t(mega_zscore_d1_selg))
mega_zscore_d3_selg <- as.data.frame(t(mega_zscore_d3_selg))

mega_zscore_d1_selg$donor <- "ET-P"

mega_zscore_d3_selg$donor <- "ET-F"

mega_zscore_selg <- rbind(mega_zscore_d1_selg,mega_zscore_d3_selg)

mega_zscore_selg <- melt(mega_zscore_selg)

mega_zscore_selg$variable <- gsub(' ','_',mega_zscore_selg$variable)

mega_zscore_selg$variable <- gsub('/','',mega_zscore_selg$variable)

mega_zscore_selg$variable <- gsub('[(]','',mega_zscore_selg$variable)

mega_zscore_selg$variable <- gsub('[)]','',mega_zscore_selg$variable)

mega_zscore_selg$variable <- factor(mega_zscore_selg$variable,levels = c('Arginine_biosynthesis','Propanoate_metabolism','Citrate_cycle_TCA_cycle','Glycolysis__Gluconeogenesis','Pyruvate_metabolism'))
ggplot(mega_zscore_selg, aes(x = variable, y = value))+ 
  
  geom_boxplot(aes(fill = donor),position=position_dodge(0.6),width=0.6,outlier.color = NA)+
  
  scale_fill_manual(values = c(cell_type_cols[1:7]))+theme_bw()+
  
  theme(panel.grid = element_blank(),
        
        panel.background = element_blank(),
        
        axis.line = element_line(),
        
        axis.text.x = element_text(size=12,color="black", face= "bold",angle=60,hjust = 1),
        
        axis.text.y = element_text(size=12, color="black", face= "bold"),
        
        plot.title = element_text(size=12,color="black", face= "bold",hjust=0))+
  
  stat_compare_means(aes(group = donor),
                     
                     label="p.signif",
                     
                     show.legend = F)

## Fig. 4I
intersect(rownames(ptim_each_zscore),rownames(copri_each_zscore)) -> id

pcopri_zscore_corr <- copri_each_zscore[c(id),]

ptim_zscore_corr <- ptim_each_zscore[c(id),]

pcopri_zscore_corr <- as.data.frame(apply(pcopri_zscore_corr,1,mean))

ptim_zscore_corr <- as.data.frame(apply(ptim_zscore_corr,1,mean))

zscore_corr <- cbind(pcopri_zscore_corr,ptim_zscore_corr)

colnames(zscore_corr) <- c('pcopri','ptim')

zscore_corr <- cbind(rownames(zscore_corr),zscore_corr)
colnames(zscore_corr)[1] <- "path_name" 


rownames(zscore_corr) <- zscore_corr$path_name

zscore_corr$anno[which(zscore_corr$anno=="metabolic pathway")] <- "Other"

b <- ggplot(zscore_corr, aes(x = ptim, y = pcopri))

b + geom_point(aes(color = anno, shape = anno))+
  
  geom_smooth(method = "lm") +
  
  #geom_rug(aes(color =outcome)) +
  
  scale_color_manual(values = c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#FDB462"))+
  
  scale_fill_manual(values = c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#FDB462"))+
  
  ggpubr::stat_cor(aes(color = anno), label.x = -1.5,label.y = c(3,2.7,2.4,2.1,1.8))+theme_bw()+
  
  theme(panel.grid=element_blank(),
        
        text=element_text(size= 15 ))



tsne_cell_cor <- cbind(mega_each_zscore_mean,copri_each_zscore_mean,ptim_each_zscore_mean,phom_each_zscore_mean)


rownames(tsne_cell_cor) <- iddd

colnames(tsne_cell_cor) <- c("mega","pcopri","ptim","phom")

crt_cor<-cor(tsne_cell_cor,method="spearman")

###pvalue
cor.mtest(tsne_cell_cor)->crt_p

cpmx<-crt_p[[1]]

rownames(cpmx)<-c("mega","pcopri","ptim","phom")

colnames(cpmx)<-c("mega","pcopri","ptim","phom")

cpmx[which(is.na(cpmx))]<-1


COLOR01<- colorRampPalette(c("#5081ff","#638dff","#b1c6ff","#c3d5ff" ,"#ffffff", "#f7cccc","#f1aeac","#eb8b8a","#e66a68"))(50)

corrplot(crt_cor,
         
         type="lower",
         #order="hclust",
         p.mat=cpmx,
         
         insig = "label_sig",
         
         #pch.col="black",
         #insig = "blank",
         sig.level = c(0.001,0.01,0.05),
         
         pch.cex = 2.5,pch.col="black",
         
         diag = FALSE,tl.srt =30,tl.col = "black",
         
         family="serif",col = COLOR01,
         
         tl.cex = 1.5,title = "",
         
         mar=c(0, 5, 2, 0),cl.cex = 0.8
)





