## Function: enrichment analysis of genes 
# -----------------------------------------
# Author: Cao Luolong, caoluolong@outlook.com
# Date created: 09-04-2021
# Date updated: 2022.07.05
# Warning: 
# @HEU.@FUDAN.
############## load packages ###########
rm(list=ls())
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.15")
# BiocManager::install("clusterProfiler")  #for enrichment
# BiocManager::install("topGO")  #for plot
# BiocManager::install("Rgraphviz")
# BiocManager::install("pathview") #for KEGG pathway
# BiocManager::install("rrvgo")
# BiocManager::install("org.Hs.eg.db") #for gene annotation.
# BiocManager::install("ReactomePA")
# devtools::install_github("jokergoo/simplifyEnrichment")

library(BiocManager)
library(clusterProfiler)
# library(R.utils)
# R.utils::setOption( "clusterProfiler.download.method",'auto')
library(topGO)
library(Rgraphviz)
# library(pathview) # library(pathviewr)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(stringr)
library(Cairo)
library(msigdbr)
library(rrvgo)
library(simplifyEnrichment)
library(GOSemSim)
library(ReactomePA)
library(DOSE)
############## load packages ###########

dir.root <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(dir.root);# kegmt<-read.gmt(paste0(dir.root,"/rawdata/msigdb.v7.5.1.symbols.gmt"))
# kegmt<-read.gmt(paste0(dir.root,"/ReactomePathways.gmt"))
# dir.root <- "F:/work_dir/AHBAenrich/enrichment_out/"
dir.work <- list.files(path = paste0(dir.root,"/enrichment_out/"),full.names = TRUE);dir.work <- dir.work[length(dir.work)]

file_name_pre <- '*total.txt$'
sig_files <- list.files(path = dir.work, pattern = file_name_pre,full.names = TRUE);out_probe_join11 <- read.table(sig_files,sep = '\t',header = T);
###########  ORA enrich  ##############
ora.enrich.all <- NULL
for(file_num in 2:ncol(out_probe_join11)){
    out_probe_join_ord <- out_probe_join11[order(-abs(out_probe_join11[,file_num]),decreasing = TRUE),c(1,file_num)]
    gene_num = round(dim(out_probe_join11)[1]/10)
    # egenes.symble <- as.character(out_probe_join11$Symbol[which(out_probe_join11$PLS3Z<=-3)])#set a threshold of gene-Z score
    egenes.symble <- as.character(out_probe_join_ord$Symbol[(dim(out_probe_join11)[1]-gene_num):dim(out_probe_join11)[1]])#set a threshhood of proporation
    file_name <- paste0('/',colnames(out_probe_join11)[file_num],'_neg')
    if (FALSE){
        # egenes.symble <- as.character(out_probe_join11$Symbol[which(out_probe_join11$PLS3Z>=1.96)])#set a threshold of gene-Z score
        egenes.symble <- as.character(out_probe_join_ord$Symbol[1:gene_num])
        file_name <- paste0('/',colnames(out_probe_join11)[file_num],'_pos')
    }# these code won't run!!!!
    egenes.entrez <- na.omit(mapIds(x = org.Hs.eg.db, keys = egenes.symble, keytype = "SYMBOL", column="ENTREZID"))
    universe.entriez <- na.omit(mapIds(x = org.Hs.eg.db, keys = out_probe_join11$Symbol, keytype = "SYMBOL", column="ENTREZID"))
    ora.gsea <- enricher(
        gene = egenes.symble,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        # universe = universe.entriez,universe = NULL,
        minGSSize = 10,
        maxGSSize = 500,
        qvalueCutoff = 0.2,
        gson = NULL,
        TERM2GENE=kegmt,
        TERM2NAME = NA)
    ora.do <- enrichDO(
        gene = egenes.entrez,
        ont = "DO",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        # universe = universe.entriez,universe = NULL,
        minGSSize = 10,
        maxGSSize = 500,
        qvalueCutoff = 0.2,
        readable = TRUE)
    ora.go <- enrichGO(gene = egenes.entrez,
         OrgDb = org.Hs.eg.db,
         keyType = "ENTREZID",
         ont = "ALL",
         pvalueCutoff = 0.01,
         pAdjustMethod = 'BH',
         # universe = universe.entriez,universe = NULL,
         qvalueCutoff = 0.2,
         minGSSize = 10,
         maxGSSize = 500,
         readable = TRUE)
    ora.wp <- enrichWP(gene = egenes.entrez,
        organism = "Homo sapiens",#get_wp_organisms()
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        # universe = universe.entriez,universe = NULL,
        minGSSize = 10,
        maxGSSize = 500,
        qvalueCutoff = 0.2)
    ora.dgn <- enrichDGN(gene = egenes.entrez,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        # universe = universe.entriez,universe = NULL,
        minGSSize = 10,
        maxGSSize = 500,
        qvalueCutoff = 0.2,
        readable = TRUE)
    ora.ncg <- enrichNCG(gene = egenes.entrez,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        # universe = universe.entriez,universe = NULL,
        minGSSize = 10,
        maxGSSize = 500,
        qvalueCutoff = 0.2,
        readable = TRUE)
    ora.kegg <- enrichKEGG(
        gene = egenes.entrez,
        organism = "hsa",
        keyType = "kegg", # one of "kegg", 'ncbi-geneid', 'ncib-proteinid' and 'uniprot'
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        pAdjustMethod = "BH",
        # universe = universe.entriez,universe = NULL,
        minGSSize = 10,
        maxGSSize = 500,
        use_internal_data = FALSE)
    ora.kegg <- setReadable(ora.kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    ora.reactome <-enrichPathway(gene = egenes.entrez,
         organism = "human",
         pvalueCutoff = 0.05,
         pAdjustMethod = "BH",
         qvalueCutoff = 0.2,
         # universe = universe.entriez,universe = NULL,
         minGSSize = 10,
         maxGSSize = 500,
         readable = T)
        ora.enrich <- list(ora.go,ora.do,ora.wp,ora.dgn,ora.ncg,ora.kegg,ora.reactome,ora.gsea);ora.enrich.all[file_num-1] <- list(ora.enrich)
    
    ###########  GO cluster  ##############    
    # simMatrix <- calculateSimMatrix(ora.go@result$ID,
    #                                 orgdb="org.Hs.eg.db",
    #                                 ont="BP",
    #                                 method="Rel")
    # scores <- setNames(-log10(ora.go@result$qvalue), ora.go@result$ID)
    # reducedTerms <- reduceSimMatrix(simMatrix,
    #                                 scores,
    #                                 threshold=0.8,
    #                                 orgdb="org.Hs.eg.db")
    # heatmapPlot(simMatrix,
    #             reducedTerms,
    #             annotateParent=TRUE,
    #             annotationLabel="parentTerm",
    #             fontsize=6)
    # treemapPlot(reducedTerms)
    # scatterPlot(simMatrix, reducedTerms)
    # 
    # mat = GO_similarity(ora.go@result$ID,ont = "BP")
    # df = simplifyGO(mat)
    # df2 <- compare_clustering_methods(mat, plot_type = "heatmap")
    ###########  GO cluster  ##############  
}
save(ora.enrich.all,list = ,file = paste0(dir.work,"/ora_enrichment1000abs.Rda"))
###########  ORA enrich  ##############

###########  GSEA enrich  ##############
gsea.enrich.all <- NULL
for(file_num in 2:ncol(out_probe_join11)){
    # file_num <- 4#here may be a loop.For file_num or id_method.
    egenes.symble <- data.frame(SYMBOL = out_probe_join11$Symbol,logFC = out_probe_join11[,file_num])#set a threshhood of gene-Z score
    gene=bitr(egenes.symble$SYMBOL,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #???в??ֻ??????ݶ?ʧ??????ENSEMBL
    gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
    gene_df <- data.frame(logFC=out_probe_join11[,file_num],SYMBOL = out_probe_join11$Symbol)
    gene_df <- merge(gene_df,gene,by="SYMBOL")
    geneList.symbol<-gene_df$logFC;names(geneList.symbol)=gene_df$SYMBOL
    geneList.entrez<-gene_df$logFC;names(geneList.entrez)=gene_df$ENTREZID
    geneList.symbol=sort(geneList.symbol,decreasing = T);geneList.entrez=sort(geneList.entrez,decreasing = T)
    gsea.gsea<-GSEA(
        geneList.symbol,
        exponent = 1,
        nPerm = 1000,
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        TERM2GENE = kegmt,
        TERM2NAME = NA,
        verbose = TRUE,
        seed = FALSE,
        by = "fgsea")
    gsea.gsea@setType <- 'MSigDb'
    gsea.do <- gseDO(
        geneList.entrez,
        exponent = 1,
        minGSSize = 10,
        maxGSSize = 500,
        eps = 1e-10,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        verbose = TRUE,
        seed = FALSE,
        by = "fgsea")
    gsea.do <- setReadable(gsea.do, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    gsea.go <- gseGO(
        geneList.symbol,
        ont = "ALL",
        OrgDb = "org.Hs.eg.db",
        keyType = "SYMBOL",
        exponent = 1,
        minGSSize = 10,
        maxGSSize = 500,
        eps = 1e-10,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        verbose = TRUE,
        seed = FALSE,
        by = "fgsea")
    gsea.wp<-gseWP(
        geneList.entrez,
        organism = "Homo sapiens",#get_wp_organisms()
        exponent = 1,
        nPerm = 1000,
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        verbose = TRUE,
        seed = FALSE,
        by = "fgsea")
    gsea.wp <- setReadable(gsea.wp, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    gsea.dgn <- gseDGN(
        geneList.entrez,
        exponent = 1,
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        verbose = TRUE,
        seed = FALSE,
        by = "fgsea")
    gsea.dgn <- setReadable(gsea.dgn, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    gsea.ncg <- gseNCG(
        geneList.entrez,
        exponent = 1,
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        verbose = TRUE,
        seed = FALSE,
        by = "fgsea")
    gsea.ncg <- setReadable(gsea.ncg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    gsea.kegg <- gseKEGG(
        geneList.entrez,
        organism = "hsa",
        keyType = "kegg",
        exponent = 1,
        minGSSize = 10,
        maxGSSize = 500,
        eps = 1e-10,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        verbose = TRUE,
        use_internal_data = FALSE,
        seed = FALSE,
        by = "fgsea")
    gsea.kegg <- setReadable(gsea.kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    gsea.reactome <- gsePathway(
        geneList.entrez,
        organism = "human",
        exponent = 1,
        minGSSize = 10,
        maxGSSize = 500,
        eps = 1e-10,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        verbose = TRUE,
        seed = FALSE,
        by = "fgsea")
    gsea.reactome2 <- setReadable(gsea.reactome, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    gsea.enrich <- list(gsea.go,gsea.do,gsea.wp,gsea.dgn,gsea.ncg,gsea.kegg,gsea.reactome,gsea.gsea);gsea.enrich.all[file_num-1] <- list(gsea.enrich)
}
save(gsea.enrich.all,list = ,file = paste0(dir.work,"/gsea_enrichment.Rda"))
###########  GSEA enrich  ##############
# ###########  visualize of diff patuways  ####
# load("E:/Project/XomicsEnrich/enrichment_out/2023_06_06_13_47PANSscz/ora_enrichment1000.Rda");gsea.enrich.all <- ora.enrich.all;scz.ora.enrich.all <- gsea.enrich.all
# load("E:/Project/XomicsEnrich/enrichment_out/2023_06_09_17_26MPasd/ora_enrichment1000.Rda");gsea.enrich.all <- ora.enrich.all;asd.ora.enrich.all <- gsea.enrich.all
enrich_type <- 4;sig_thresh <- 5e-2#;gsea.enrich.all <- ora.enrich.all
data1 <- gsea.enrich.all[[1]][[enrich_type]]@result %>%filter(qvalue<=sig_thresh);data2 <- gsea.enrich.all[[2]][[enrich_type]]@result %>%filter(qvalue<=sig_thresh)
test<- full_join(select(data1,ID,qvalue,Description),
                 select(data2,ID,qvalue,Description),by = "ID",suffix = c("", ".2"))
test$Description[is.na(test$Description)] <- test$Description.2[is.na(test$Description)];test<- select(test,-ncol(test))

for(file_num in 3:length(gsea.enrich.all)){
    # for(file_num in 9:14){
    data2 <- gsea.enrich.all[[file_num]][[enrich_type]]@result %>%filter(qvalue<=sig_thresh)
    test<- full_join(test,
                     select(data2,ID,qvalue,Description),by = "ID",suffix = c("", ".2"))
    test$Description[is.na(test$Description)] <- test$Description.2[is.na(test$Description)];test<- select(test,-ncol(test))
}
test <- subset(test,grepl("depressi|bipol|schiz|autis|alcoh",Description))
# SCZ
test <- subset(test,grepl("dopamin|glutamat|neurodevelop|gaba|gamma|immune|cytoskel|synaptic",Description))
# test <- subset(test,grepl("dopamin|glutamat|neurodevelop|gaba|gamma|immune|inflamma|epigenetic|oxidative",Description))
# # ASD
# test <- subset(test,grepl("synaptic|neural|epigenetic|immune|metaboli|mitochondri",Description))

test[is.na(test)] <- 1;colnames(test)[4:16] <- colnames(out_probe_join11)[3:ncol(out_probe_join11)]
test.new <- select(test,-3,-1)
rownames(test.new) <- test$Description;colnames(test.new) <- colnames(out_probe_join11)[2:ncol(out_probe_join11)]
test.new <- -log10(test.new)
ComplexHeatmap::pheatmap(test.new,color = colorRampPalette(c("white","orange"))(50),
                         cluster_rows = TRUE,cluster_cols = F,main = gsea.enrich.all[[1]][[enrich_type]]@ontology)#ontology,setType

library(wordcloud)
wordcloud(words = test$Description,freq = -log10(test$max_Z_PPI),max.words = 100,col=brewer.pal(8,"Dark2"))
library(wordcloud2)
data1 <- data.frame(test$Description,-log10(test$max_Z_PPI));data1 <- data1[order(data1$X.log10.test.max_Z_PPI.,decreasing = TRUE),];data1 <- data1[1:1000,]
data2 <- data.frame(test$Description,-log10(test$PLSg));data2 <- data2[order(data2$X.log10.test.PLSg.,decreasing = TRUE),];data2 <- data2[1:1000,]
# wordcloud(words = data1$test.Description,freq = data1$X.log10.test.max_Z_PPI.,min.freq=0,col=brewer.pal(8,"Dark2"))
wordcloud2(data = data1, shuffle = F,color = brewer.pal(8,"Dark2"))
wordcloud2(data = data2, shuffle = F,color = brewer.pal(8,"Dark2"))


#####   heatmap
library("pheatmap")
X40corr <- read.table("clipboard",header=T)
rownames(X40corr) <- colnames(X40corr)
pheatmap(X40corr,color = colorRampPalette(c("blue","white","orange"))(paletteLength),breaks = c(seq(-1,1,length.out= paletteLength)),
         fontsize = 10,angle_col = 45,cluster_rows = TRUE,cluster_cols = T,main = "Correlation among PET maps")#ontology,setType
library(showtext)
windowsFonts();font_families();
#####  line chart
library(reshape2)#融合数据
library(ggplot2)#绘图工具
require(cowplot)
data = read.table("clipboard",header = T);data <- as.data.frame(t(data[3:7,]));data$method <- rownames(data);data <- melt(data,id = "method")
colnames(data) <- c("method","top_set","genes")#更改列名
l2 <- ggplot(data = data,aes(x=method,y=genes,group = top_set,color=top_set,shape=top_set))+
    geom_point(size=3)+
    geom_line(linetype=1,size=1)+
    xlab("method")+#横坐标名称
    ylab("# of overlap ASD genes")+#纵坐标名称
    theme_bw() +#去掉背景灰色
    theme(panel.grid.major=element_line(colour=NA),
          panel.grid.minor = element_blank(),#以上theme中代码用于去除网格线且保留坐标轴边框
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          axis.title = element_text(size = 15),#设置坐标轴名称
          axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75,size = 10),#横坐标倾斜
          axis.text.y = element_text(angle = 0,vjust = 0.85,hjust = 0.75,size = 10),#纵坐标倾斜
          text = element_text(family = "STXihei"),#设置中文字体的显示
          legend.position = c(.035,.73),#更改图例的位置，放至图内部的左上角
          legend.box.background = element_rect(color="black"))#为图例添加边框线
    #+scale_x_continuous(limits = c(2000,2018),breaks = seq(2000,2018,1))更改横坐标刻度值
plot_grid(l1,l2,labels = LETTERS[1:2], nrow = 2)

#########  scattter plot
sig_files <- "E:/Project/XomicsEnrich/enrichment_out/2023_06_06_13_47PANSscz/map2whole308_PLS_boot_total.txt" # SCZ
out_probe_join11 <- read.csv(sig_files,sep = '\t')#here may be a loop.
test <- out_probe_join11[,c(7,2,10,13,1)];test[,2] <- 0.5*test[,2]/max(abs(test[,2]));test[,3] <- 0.5*test[,3]/max(abs(test[,3]));test[,4] <- 0.5*test[,4]/max(abs(test[,4]))
"#FB8072""#80B1D3"
# sig_Z <- 1.1;test$color <- "else";test$color[intersect(which(test$PLSboot_gene/test$R>sig_Z),which(abs(test$PLSboot_gene)>0.4))] <- "high";test$color[intersect(which(test$R/test$PLSboot_gene>sig_Z),which(abs(test$R)>0.4))] <- "low"
# df1 <- test %>%filter(color == "low");df1 <-df1[order(abs(df1$R),decreasing = T)[1:20],]
# df2 <- test %>%filter(color == "high");df2 <-df2[order(abs(df2$PLSboot_gene),decreasing = T)[1:20],]
sig_Z <- 1000;test$color <- "else";test$diff_ord <- abs(rank(test$PLSboot_gene)-rank(test$R))#从小到大
test$color[intersect(which(test$diff_ord>sig_Z),order(abs(test$PLSboot_gene),decreasing = T)[1:sig_Z])] <- "high"
test$color[intersect(which(test$diff_ord>sig_Z),order(abs(test$R),decreasing = T)[1:sig_Z])] <- "low"
df1 <- test %>%filter(color == "high");df1 <-df1[order(abs(df1$PLSboot_gene),decreasing = T)[1:20],]
df2 <- test %>%filter(color == "low");df2 <-df2[order(abs(df2$R),decreasing = T)[1:20],]
p1 <- ggplot(test,aes(x = PLSboot_gene,y=R))+
    geom_point(aes(colour = color,size = 0.1),alpha=0.9)+
    scale_colour_manual(values=c("low" = "blue","else" = "grey","high" = "orange"))+
    theme_bw()+
    # labs(title = "abc")+
    scale_x_continuous(limits = c(-0.5, 0.5))+
    theme(panel.grid.major=element_line(colour=NA),
          axis.text = element_text(angle = 0,vjust = 0.85,hjust = 0.75,size = 10))+# family = "serif",
    geom_abline(intercept = 0, slope = 1)+
    ggrepel::geom_text_repel(aes(label=Symbol),rbind(df1,df2),max.overlaps =50)
# test$color <- "else";test$color[intersect(which(test$max_Z_PPI/test$R>sig_Z),which(abs(test$max_Z_PPI)>0.4))] <- "high";test$color[intersect(which(test$R/test$max_Z_PPI>sig_Z),which(abs(test$R)>0.4))] <- "low"
# df1 <- test %>%filter(color == "low");df1 <-df1[order(abs(df1$R),decreasing = T)[1:20],]
# df2 <- test %>%filter(color == "high");df2 <-df2[order(abs(df2$max_Z_PPI),decreasing = T)[1:20],]
test$color <- "else";test$diff_ord <- abs(rank(test$max_Z_PPI)-rank(test$R))#从小到大
test$color[intersect(which(test$diff_ord>sig_Z),order(abs(test$max_Z_PPI),decreasing = T)[1:sig_Z])] <- "high"
test$color[intersect(which(test$diff_ord>sig_Z),order(abs(test$R),decreasing = T)[1:sig_Z])] <- "low"
df1 <- test %>%filter(color == "high");df1 <-df1[order(abs(df1$max_Z_PPI),decreasing = T)[1:20],]
df2 <- test %>%filter(color == "low");df2 <-df2[order(abs(df2$R),decreasing = T)[1:20],]
p2 <- ggplot(test,aes(x = max_Z_PPI,y=R))+
    geom_point(aes(colour = color,size = 0.1),alpha=0.9)+
    scale_colour_manual(values=c("low" = "blue","else" = "grey","high" = "orange"))+
    theme_bw()+
    # labs(title = "abc")+
    scale_x_continuous(limits = c(-0.5, 0.5))+
    theme(panel.grid.major=element_line(colour=NA),
          axis.text = element_text(angle = 0,vjust = 0.85,hjust = 0.75,size = 10))+# family = "serif",
    geom_abline(intercept = 0, slope = 1)+ 
    ggrepel::geom_text_repel(aes(label=Symbol),rbind(df1,df2),max.overlaps =50)
# test$color <- "else";test$color[intersect(which(test$PLSg/test$R>sig_Z),which(abs(test$PLSg)>0.4))] <- "high";test$color[intersect(which(test$R/test$PLSg>sig_Z),which(abs(test$R)>0.4))] <- "low"
# df1 <- test %>%filter(color == "low");df1 <-df1[order(abs(df1$R),decreasing = T)[1:20],]
# df2 <- test %>%filter(color == "high");df2 <-df2[order(abs(df2$PLSg),decreasing = T)[1:20],]
test$color <- "else";test$diff_ord <- abs(rank(test$PLSg)-rank(test$R))#从小到大
test$color[intersect(which(test$diff_ord>sig_Z),order(abs(test$PLSg),decreasing = T)[1:sig_Z])] <- "high"
test$color[intersect(which(test$diff_ord>sig_Z),order(abs(test$R),decreasing = T)[1:sig_Z])] <- "low"
df1 <- test %>%filter(color == "high");df1 <-df1[order(abs(df1$PLSg),decreasing = T)[1:20],]
df2 <- test %>%filter(color == "low");df2 <-df2[order(abs(df2$R),decreasing = T)[1:20],]
p3 <- ggplot(test,aes(x = PLSg,y=R))+
    geom_point(aes(colour = color,size = 0.1),alpha=0.9)+
    scale_colour_manual(values=c("low" = "blue","else" = "grey","high" = "orange"))+
    theme_bw()+
    # labs(title = "abc")+
    scale_x_continuous(limits = c(-0.5, 0.5))+
    theme(panel.grid.major=element_line(colour=NA),
          axis.text = element_text(angle = 0,vjust = 0.85,hjust = 0.75,size = 10))+# family = "serif",
    geom_abline(intercept = 0, slope = 1)+ 
    ggrepel::geom_text_repel(aes(label=Symbol),rbind(df1,df2),max.overlaps =50)

plot_grid(p1,p2,p3,labels = LETTERS[1:3], nrow = 1)
# add_el <- theme_grey() + theme(text = element_text(family = "Times"))
