setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

pcg <- c("plyr","stringr","magrittr","RColorBrewer","ggplot2","Mfuzz","readxl","pheatmap","dplyr","ggpubr","ggbiplot","ComplexHeatmap","reshape2","openxlsx","circlize","ggrepel","SeqKnn","limma","umap")
sapply(pcg,library,character.only=T)

`%notin%` <- Negate(`%in%`)
volcano.plot <- function(data,#data.frame
                             group1,
                             group2,
                             label1,
                             label2,
                             paired = F,
                             pvalue= 0.05,
                             P.adjust= F,
                             var.equal = F,
                             select_gene = NULL,
                             fold_change1 = 1.2,
                             fold_change2 = 1.2,
                             label_size = 3){
  df1 <- data
  #填充NA
  df1[is.na(df1)] <- 0
  # df1[is.na(df1)]=min(df1,na.rm = T)
  #计算foldchange
  df1$`log2(foldchange)` <- apply(df1,1, function(x) log2((mean(x[group1],na.rm = T)/mean(x[group2],na.rm = T))))
  
  if(paired){
    P_value <- apply(df1,
                     1,
                     function(y) {
                       
                       p_try = tryCatch(t.test(y[group1],
                                               y[group2],
                                               paired = paired,
                                               var.equal = var.equal)$p.value,
                                        error = function(x) NA)
                     })
  }else{
    P_value <- apply(df1,
                     1,
                     function(y) {
                       
                       p_try = tryCatch(t.test(y[group1],
                                               y[group2],
                                               paired = F,
                                               var.equal = var.equal)$p.value,
                                        error = function(x) NA)
                     })  
  }
  #给data增加p-value列
  df1$P_value<- P_value
  df1$P_value_adjust<-p.adjust(df1$P_value, method="BH")
  
  
  if(P.adjust){
    df1$threshold = factor(ifelse(df1$P_value_adjust < pvalue & (df1$`log2(foldchange)` > log2(fold_change1) | df1$`log2(foldchange)` < -log2(fold_change2)),ifelse(df1$`log2(foldchange)` > log2(fold_change1),'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
    
    df1$pro <-  rownames(df1)
    up <- subset(df1, df1$P_value_adjust < pvalue & df1$`log2(foldchange)` > log2(fold_change1))
    down <- subset(df1, df1$P_value_adjust < pvalue & df1$`log2(foldchange)` < -log2(fold_change2))
    print(paste0("down:",nrow(down)))
    print(paste0("up:",nrow(up)))
    write.csv(up,file = paste0(label1,"_",label2, "_up_volcano.csv"))
    write.csv(down,file = paste0(label1,"_",label2, "_dw_volcano.csv"))
    write.csv(df1,file = paste0(label1,"_",label2, "_all_volcano.csv"))
    differprot <- list(up=row.names(up),down=row.names(down))
    
    if(!is.null(select_gene)) {data_text = df1[select_gene,]}
    else {data_text = subset(df1, P_value_adjust < pvalue & (df1$`log2(foldchange)` > log2(fold_change1) | df1$`log2(foldchange)` < -log2(fold_change2)))}
    
    pp <- ggplot(df1,aes(x=`log2(foldchange)`,y=-log10(P_value_adjust),color=threshold))+
      geom_point()+
      scale_color_manual(values=c('Up'="#DC143C",'Down'="#0A9731",'NoSignifi'="#808080"))+#确定点的颜色
      geom_text_repel(
        # data = subset(df1, P_value_adjust < pvalue & (df1$`log2(foldchange)` > log2(fold_change1) | df1$`log2(foldchange)` < -log2(fold_change2))),
        data = data_text,
        aes(label = pro),
        colour = "black",  ###标注出来的基因的颜色
        size = label_size,
        max.iter=2500,
        segment.color = "black", show.legend = FALSE )+
      theme_bw()+#修改图片背景
      theme(
        legend.title = element_blank()#不显示图例标题
      )+
      ylab("-log10 (adjust P value)")+#修改y轴名称
      xlab(paste("log2 (fold change)"))+#修改x轴名称
      geom_vline(xintercept=c(-log2(fold_change2),log2(fold_change1)),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
      geom_hline(yintercept = -log10(pvalue),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05
    #dev.off()
    ggsave(paste0(label1,"_",label2, "_volcano_id.pdf"),plot =pp ,width=8,height=6,device = NULL)
  }else{
    df1$threshold = factor(ifelse(df1$P_value < pvalue & (df1$`log2(foldchange)` > log2(fold_change1) | df1$`log2(foldchange)` < -log2(fold_change2)),ifelse(df1$`log2(foldchange)` > log2(fold_change1),'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
    df1$pro <-  rownames(df1)
    #画图
    if(!is.null(select_gene)) {data_text = df1[select_gene,]}
    else {data_text = subset(df1, P_value < pvalue & (df1$`log2(foldchange)` > log2(fold_change1) | df1$`log2(foldchange)` < -log2(fold_change2)))}
    
    pp <- ggplot(df1,aes(x=`log2(foldchange)`,y=-log10(P_value),color=threshold))+
      
      geom_point()+
      scale_color_manual(values=c('Up'="#DC143C",'Down'="#0A9731",'NoSignifi'="#808080"))+#确定点的颜色
      geom_text_repel(
        data = data_text,
        # data = subset(df1, P_value < pvalue & (df1$`log2(foldchange)` > log2(fold_change1) | df1$`log2(foldchange)` < -log2(fold_change2))),
        aes(label = pro),
        colour = "black", ###标注出来的基因的颜色
        size = label_size,
        max.iter=3000,
        segment.color = "black", show.legend = FALSE )+
      
      theme_bw()+#修改图片背景
      theme(
        legend.title = element_blank()#不显示图例标题
      )+
      ylab("-log10 (P value)")+#修改y轴名称
      xlab(paste("log2 (fold change)"))+#修改x轴名称
      geom_vline(xintercept=c(-log2(fold_change2),log2(fold_change1)),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
      geom_hline(yintercept = -log10(pvalue),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05
    ggsave(paste0(label1,"_",label2, "_volcano_id.pdf"),plot =pp ,width=12,height=8,device = NULL)
    
    #存储上调和下调蛋白信息
    up <- subset(df1, df1$P_value < pvalue & df1$`log2(foldchange)` > log2(fold_change1))
    down <- subset(df1, df1$P_value < pvalue & df1$`log2(foldchange)` < -log2(fold_change2))
    print(paste0("down:",nrow(down)))
    print(paste0("up:",nrow(up)))
    write.csv(up,file = paste0(label1,"_",label2, "_up_volcano.csv"))
    write.csv(down,file = paste0(label1,"_",label2, "_dw_volcano.csv"))
    write.csv(df1,file = paste0(label1,"_",label2, "_all_volcano.csv"))
    differprot <- list(up=row.names(up),down=row.names(down))
    
  }
  return(differprot)
}

############ PCA
plot.pca <- function(data,type,title="",ptColors=NULL,label2=NULL,width=12,height=8,ellise_type=NULL){ 
  M <- t(data)
  M[is.na(M)] <- 0
  M <- apply(M,2,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})
  clnames <- row.names(data)
  set.seed(10)
  m1 <- prcomp(M);
  Y  <- scale(M, m1$center, m1$scale) %*% m1$rotation 
  Y  <- Y[,c(1,2)]
  
  Y <- data.frame(Y,type);
  colnames(Y) <- c("PC1","PC2","label")
  eigs <- m1$sdev^2
  percentages <- eigs[1:2] / sum(eigs)
  
  p <- ggplot(Y, aes(x=PC1, y=PC2, colour=label)) + geom_point(size=4)
  p <- p + theme(  panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.text = element_text(size = 30,color = "black"),
                   panel.border = element_blank(),
                   axis.line.x = element_line(color="black", size = 0.25),
                   axis.line.y = element_line(color="black", size = 0.25),
                   plot.title   = element_text(size=30),
                   axis.title   =element_text(size=30),
                   panel.background = element_blank())
  
  strLabx <- sprintf("PC1(%4.2f%%)",percentages[1]*100)
  p <- p +  labs(x =strLabx,y = sprintf("PC2(%4.2f%%)",percentages[2]*100),
                 title =sprintf("PCA:%d features",length(clnames)))
  if(!is.null(ptColors)){
    p <- p +   scale_colour_manual(values=ptColors)
  }
  if(!is.null(label2)){
    p <- p +   geom_text(aes(label=label2,vjust = -0.8, hjust = 0.5,size=0.5),show.legend = FALSE)
  }
  if(!is.null(ellise_type)){
    p <- p +stat_ellipse(level = 0.95, show.legend = F)
  }
  ggsave(paste0(title,"_PCA.pdf"),plot =p ,width=width,height=height,device = NULL)
}

get.color <- function(n){
  if(n<10){
    set.seed(21)
    m <- sample(9,n)
    color <- brewer.pal(9,"Set1")[m]
  }else{
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    set.seed(10)
    color=sample(col_vector, n)
  }
  return(color)
}

paired_col <- c(brewer.pal(11,"Paired"))
TNM_col <- c(brewer.pal(11,"RdYlBu"))
stage_col <- c(brewer.pal(11,"PRGn"))


########  TMT 2022.2.26 血浆血清 batch7
batch7_info <- read.xlsx("./data/TMT_batch7_sample_info_20220226-1238.xlsx",sheet = "Sheet1")
raw_batch7 <-read_xlsx("./data/b7.xlsx")     # 后搜的第7个batch的数据
raw7 <- raw_batch7[raw_batch7$`Protein FDR Confidence: Combined`=="High",]

#### 第7个batch的数据
df_b7<-raw7[,grepl("Abundance",names(raw7))]
df_b7<-as.data.frame(df_b7)
row.names(df_b7)<-paste0(raw7$Accession,"^",raw7$`Gene Symbol`)
df_b7$protein <- rownames(df_b7)

colnames(df_b7) <- gsub("Ratio: \\(","Ratio: (F7, ",colnames(df_b7))
colnames(df_b7) <- gsub("126","F7, 126",colnames(df_b7))
colnames(df_b7) <- gsub("\\(Grouped\\): ","(Grouped): F7, ",colnames(df_b7))
colnames(df_b7) <- gsub("F7, F7, ","F7, ",colnames(df_b7))

con_id_b7 <- grep("CON.*",rownames(df_b7),value=TRUE)
df2_b7 <- df_b7[rownames(df_b7)[!rownames(df_b7) %in% con_id_b7], ]
df2_b7<-df2_b7[apply(df2_b7,1,function(x){sum(!is.na(x))})>0,]

df3_b7 <- df2_b7[,grep("Ratio|Grouped...F[0-9]..126",names(df2_b7))]

name=str_extract(names(df3_b7),"F[0-9]..[0-9A-Z]*")
names(df3_b7)= gsub("\\, ","_",name)

df4_b7 <- df3_b7[,-grep("126",names(df3_b7))]
df4_b7 <- df4_b7[apply(df4_b7, 1, function(x) {sum(!is.na(x))>0}),]


### 离群值处理
df4_b7_iqr <- df4_b7
df4_b7_iqr.matrix <- as.matrix(df4_b7_iqr)
outline <- (fivenum(df4_b7_iqr.matrix)[4]-fivenum(df4_b7_iqr.matrix)[2])*2+fivenum(df4_b7_iqr.matrix)[4]
df4_b7_iqr[df4_b7_iqr>outline] <- outline

group_plasma <- batch7_info[batch7_info$Group.b=="plasma",]$MS_ID
group_serum <- batch7_info[batch7_info$Group.b=="serum",]$MS_ID

df2_Gps <- df4_b7_iqr[,c(group_plasma,group_serum)]
df2_Gps <- as.data.frame(df2_Gps)
df2_Gps_diff_protein <- volcano.plot(df2_Gps,c(1:length(group_plasma)),c((length(group_plasma)+1):dim(df2_Gps)[2]),"plasma","serum",P.adjust = F,fold_change1 =  1.2,fold_change2 = 1.2)

##### 2022.2.26 proteinGroups
raw_proteinGroup <- read.csv("./data/proteinGroups(1).txt",sep = "\t",check.names = FALSE)
colname_1 <- grep("LFQ intensity ",colnames(raw_proteinGroup),value=TRUE)
colname_2 <- colname_1[grep("Plasma|Serum",colname_1)]

df1_proteinGroup <- raw_proteinGroup[,c("Majority protein IDs","Gene names",colname_2)]

rownames(df1_proteinGroup) <- paste0(df1_proteinGroup$`Majority protein IDs`,"^",df1_proteinGroup$`Gene names`)

df2_proteinGroup <- df1_proteinGroup[,c(3:dim(df1_proteinGroup)[2])]

con_id_proteing <- grep("CON.*",rownames(df2_proteinGroup),value=TRUE)
df3_proteinGroup <- df2_proteinGroup[rownames(df2_proteinGroup)[!rownames(df2_proteinGroup) %in% con_id_proteing], ]
rev_id_proteing <- grep("REV_.*",rownames(df3_proteinGroup),value=TRUE)
df3_proteinGroup <- df3_proteinGroup[rownames(df3_proteinGroup)[!rownames(df3_proteinGroup) %in% rev_id_proteing], ]

df3_proteinGroup <- df3_proteinGroup[apply(df3_proteinGroup,1,function(x){sum(!is.na(x))})>0,]


group_plasma_1 <- grep("Plasma",colnames(df2_proteinGroup),value=TRUE)
group_serum_1 <- grep("Serum",colnames(df2_proteinGroup),value=TRUE)

df2_proteinGroups_ps <- df3_proteinGroup[,c(group_plasma_1,group_serum_1)]
df2_proteinGroups_ps <- as.data.frame(df2_proteinGroups_ps)
df2_Gps_pG_diff_protein <- volcano.plot(df2_proteinGroups_ps,c(1:length(group_plasma_1)),c((length(group_plasma_1)+1):dim(df2_proteinGroups_ps)[2]),"proteinGroup_plasma","serum",P.adjust = F,fold_change1 =  1.2,fold_change2 = 1.2)

################# 新冠数据分析 start
sample_info <- read.xlsx("./data/sTable1_20220217-2046_shiyq_edit.xlsx",sheet = "Sheet2")
raw <-read_xlsx("./data/o_1-(1).xlsx")  # 先搜的第一个batch的数据
raw1 <- raw[raw$`Protein FDR Confidence: Combined`=="High",]

raw_1 <-read_xlsx("./data/b2.xlsx")     # 后搜的第二个batch的数据
raw2 <- raw_1[raw_1$`Protein FDR Confidence: Combined`=="High",]

raw_3 <-read_xlsx("./data/b3.xlsx")     # 后搜的第三个batch的数据
raw3 <- raw_3[raw_3$`Protein FDR Confidence: Combined`=="High",]

raw_4 <-read_xlsx("./data/b4.xlsx")     # 后搜的第三个batch的数据
raw4 <- raw_4[raw_4$`Protein FDR Confidence: Combined`=="High",]

raw_5 <-read_xlsx("./data/CAB20220222sunr_CVDO_TMT16plex_60min_b5_.xlsx")     # 后搜的第三个batch的数据
raw5 <- raw_5[raw_5$`Protein FDR Confidence: Combined`=="High",]

raw_6 <-read_xlsx("./data/CAA20220222sunr_CVDO_TMT16plex_60min_b6_.xlsx")     # 后搜的第三个batch的数据
raw6 <- raw_6[raw_6$`Protein FDR Confidence: Combined`=="High",]

#### 第一个batch的数据
df1<-raw1[,grepl("Abundance",names(raw1))]
df1<-as.data.frame(df1)
row.names(df1)<-paste0(raw1$Accession,"^",raw1$`Gene Symbol`)
df1$protein <- rownames(df1)
colnames(df1) <- gsub("F2, ","F1, ",colnames(df1))
#### 第二个batch的数据
d11<-raw2[,grepl("Abundance",names(raw2))]
d11<-as.data.frame(d11)
row.names(d11)<-paste0(raw2$Accession,"^",raw2$`Gene Symbol`)
d11$protein <- rownames(d11)

colnames(d11) <- gsub("Ratio: \\(","Ratio: (F2, ",colnames(d11))
colnames(d11) <- gsub("126","F2, 126",colnames(d11))
colnames(d11) <- gsub("\\(Grouped\\): ","(Grouped): F2, ",colnames(d11))
colnames(d11) <- gsub("F2, F2, ","F2, ",colnames(d11))
#### 第三个batch的数据
d22<-raw3[,grepl("Abundance",names(raw3))]
d22<-as.data.frame(d22)
row.names(d22)<-paste0(raw3$Accession,"^",raw3$`Gene Symbol`)
d22$protein <- rownames(d22)

colnames(d22) <- gsub("Ratio: \\(","Ratio: (F3, ",colnames(d22))
colnames(d22) <- gsub("\\(Grouped\\): ","(Grouped): F3, ",colnames(d22))
colnames(d22) <- gsub("126","F3, 126",colnames(d22))
colnames(d22) <- gsub("F3, F3, ","F3, ",colnames(d22))
#### 第四个batch的数据
d33<-raw4[,grepl("Abundance",names(raw4))]
d33<-as.data.frame(d33)
row.names(d33)<-paste0(raw4$Accession,"^",raw4$`Gene Symbol`)
d33$protein <- rownames(d33)

colnames(d33) <- gsub("Ratio: \\(","Ratio: (F4, ",colnames(d33))
colnames(d33) <- gsub("\\(Grouped\\): ","(Grouped): F4, ",colnames(d33))
colnames(d33) <- gsub("126","F4, 126",colnames(d33))
colnames(d33) <- gsub("F4, F4, ","F4, ",colnames(d33))
#### 第五个batch的数据
d44<-raw5[,grepl("Abundance",names(raw5))]
d44<-as.data.frame(d44)
row.names(d44)<-paste0(raw5$Accession,"^",raw5$`Gene Symbol`)
d44$protein <- rownames(d44)

colnames(d44) <- gsub("Ratio: \\(","Ratio: (F5, ",colnames(d44))
colnames(d44) <- gsub("\\(Grouped\\): ","(Grouped): F5, ",colnames(d44))
colnames(d44) <- gsub("126","F5, 126",colnames(d44))
colnames(d44) <- gsub("F5, F5, ","F5, ",colnames(d44))
#### 第六个batch的数据
d55<-raw6[,grepl("Abundance",names(raw6))]
d55<-as.data.frame(d55)
row.names(d55)<-paste0(raw6$Accession,"^",raw6$`Gene Symbol`)
d55$protein <- rownames(d55)

colnames(d55) <- gsub("Ratio: \\(","Ratio: (F6, ",colnames(d55))
colnames(d55) <- gsub("\\(Grouped\\): ","(Grouped): F6, ",colnames(d55))
colnames(d55) <- gsub("126","F6, 126",colnames(d55))
colnames(d55) <- gsub("F6, F6, ","F6, ",colnames(d55))

# 数据合并
df_22 <- merge(df1,d11,by = "protein",all = TRUE)
df_22 <- merge(df_22,d22,by = "protein",all = TRUE)
df_22 <- merge(df_22,d33,by = "protein",all = TRUE)
df_22 <- merge(df_22,d44,by = "protein",all = TRUE)
df_22 <- merge(df_22,d55,by = "protein",all = TRUE)

rownames(df_22) <- df_22$protein
df_22 <- df_22[,c(2:dim(df_22)[2])]
con_id <- grep("CON.*",rownames(df_22),value=TRUE)
new_df2 <- df_22[rownames(df_22)[!rownames(df_22) %in% con_id], ]
new_df2<-new_df2[apply(new_df2,1,function(x){sum(!is.na(x))})>0,]

df3 <- new_df2[,grep("Ratio|Grouped...F[0-9]..126",names(new_df2))]

name=str_extract(names(df3),"F[0-9]..[0-9A-Z]*")
names(df3)= gsub("\\, ","_",name)

df4 <- df3[,-grep("126",names(df3))]
pool<- df3[,grep("126",names(df3))]

pool <- pool[apply(pool, 1, function(x) {sum(!is.na(x))>0}),]
df4 <- df4[apply(df4, 1, function(x) {sum(!is.na(x))>0}),]

################# pool样本CV   FigS1 C
tmp.cv <- apply(pool, 1 , function(x){sd(x,na.rm = T)/mean(x,na.rm =T)})
tmp.cv <-as.data.frame(tmp.cv)
df_pool.cv <- data.frame(cv=tmp.cv,sample=rep("pool",each=nrow(tmp.cv)))
colnames(df_pool.cv) <- c("cv","sample")
#删除cv为NA的行
df_pool.cv <- df_pool.cv[!is.na(df_pool.cv$cv),]
# write.csv(df_pool.cv,file = "pool_CV.csv")
###########plot violin
n_fun <- function(x){
  return(data.frame(y = median(x), label = paste0("median =",round(median(x),3))))
}
p<-ggplot(df_pool.cv, aes(x = sample, y=cv,color=sample)) +
  xlab("")+
  geom_violin(trim=FALSE)+
  scale_fill_manual(values=brewer.pal(12,"Set3")[c(1:10)])+
  geom_boxplot(width=0.1)+
  theme(legend.direction = 'horizontal',legend.position = 'top',panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  theme(axis.text = element_text(size = 15,colour = "black"),text = element_text(size = 15,colour = "black"))+
  theme(axis.text.x = element_text( hjust = 1,angle = 45))+
  theme(legend.position = "none")+
  stat_summary(fun.data =n_fun, geom="text")
ggsave("FigS1 C.pdf",plot =p ,device = NULL,width = 4, height = 4)

new_df4<-df4[apply(df4,1,function(x){sum(!is.na(x))})>0.2*dim(df4)[2],]

df4_seqKNN <- SeqKNN(new_df4,10)

### 离群值处理
# df4_iqr <- new_df4
df4_iqr <- df4_seqKNN
df4_iqr.matrix <- as.matrix(df4_iqr)
outline <- (fivenum(df4_iqr.matrix)[4]-fivenum(df4_iqr.matrix)[2])*2+fivenum(df4_iqr.matrix)[4]
df4_iqr[df4_iqr>outline] <- outline

############### 删除血清 血浆里面的差异蛋白
del_plasma_protein <- unique(c(sapply(sapply(df2_Gps_pG_diff_protein$up, function(x) str_split(x,"\\^")[[1]][1]),function(x) str_split(x,";")[[1]][1]),
                               sapply(sapply(df2_Gps_pG_diff_protein$down, function(x) str_split(x,"\\^")[[1]][1]),function(x) str_split(x,";")[[1]][1]),
                               sapply(sapply(df2_Gps_diff_protein$up, function(x) str_split(x,"\\^")[[1]][1]),function(x) str_split(x,";")[[1]][1]),
                               sapply(sapply(df2_Gps_diff_protein$down, function(x) str_split(x,"\\^")[[1]][1]),function(x) str_split(x,";")[[1]][1])))

new_df4_iqr <- df4_iqr
rownames(new_df4_iqr) <- sapply(rownames(new_df4_iqr), function(x) str_split(x,"\\^")[[1]][1])
df4_plasma_overlap_protein <- intersect(rownames(new_df4_iqr),del_plasma_protein)
df4_plasma_TMT_protein <- new_df4_iqr[rownames(new_df4_iqr) %notin% df4_plasma_overlap_protein,]
###### 2022.2.28 Gene Symbol
genesymbol_0228 <- read.csv("./data/gn751.txt",sep = "\t")

rownames(df4_plasma_TMT_protein) <- paste0(rownames(df4_plasma_TMT_protein),"_",genesymbol_0228$To[match(rownames(df4_plasma_TMT_protein),genesymbol_0228$From)])

new_tmt_sample_info <- read.xlsx("./data/TMT_sample_info_20220225.xlsx")

##### 技术重复 FigS1 B
df_f.name <- new_tmt_sample_info$SampleID[match(colnames(df4_plasma_TMT_protein),new_tmt_sample_info$MS_ID)]
repB <- df_f.name[grepl("_rep",df_f.name)]
repB1 <- unique(gsub("_rep","",repB))

pdf("FigS1 B.pdf",height =8,width = 8)
par(mfrow = c(2, 2))
for (nam in repB1) {
  rep_sample <- df_f.name[grepl(nam,df_f.name)]
  tmp_info <- new_tmt_sample_info[new_tmt_sample_info$SampleID %in% rep_sample,]
  
  df.r <- df4_plasma_TMT_protein[,tmp_info$MS_ID]
  
  cor1 <- as.numeric(df.r[,1]) %>% scale()
  cor2 <- as.numeric(df.r[,2]) %>% scale()
  r <- cor(cor1, cor2, use = "pairwise.complete.obs",method = "pearson")   
  smoothScatter(cor1, cor2, nrpoints = 100,cex = 2,
                colramp = colorRampPalette(c(blues9,"orange", "red")),
                main = nam, xlab = "repA", ylab = "repB")
  abline(lm(cor1 ~ cor2), col="red", lwd=2, lty=2)
  text(-max(cor1,na.rm = T)*0.2,max(cor2,na.rm = T)*0.8,labels =paste0( "r =", as.character(round(r,4))),cex = 2)
  
}
dev.off()

data_diff <- df4_plasma_TMT_protein

########## 差异分析 分组
group5 <- new_tmt_sample_info[new_tmt_sample_info$Group.b=="5",]$MS_ID
group6 <- new_tmt_sample_info[new_tmt_sample_info$Group.b=="6",]$MS_ID

df2_G56 <- data_diff[,c(group5,group6)]
df2_G56 <- as.data.frame(df2_G56)
df2_G56_diff_protein_pvalue <- volcano.plot(df2_G56,c(1:length(group5)),c((length(group5)+1):dim(df2_G56)[2]),"Group_5","6_pvalue",P.adjust = F,fold_change1 =  1.2,fold_change2 = 1.2)


###### 2022.2.28 line chart MainFig C
# 对每个蛋白的中位数取值做成一条线, 然后把在group 56里面的差异蛋白(一共是167个蛋白)的表达情况全都做在一张图上(横坐标是group 5, 6, 1的顺序, 纵坐标是Normalized protein ratio), 类似于mfuzz的那种图. 或者你如果不会做, 看看能不能导入mfuzz, 直接做一下. 应该是如下效果
# df4_plasma_TMT_protein
# df2_G56_diff_protein_pvalue
df4_G45_diff_protein <- df4_plasma_TMT_protein[c(df2_G56_diff_protein_pvalue$up,df2_G56_diff_protein_pvalue$down),]

temp_heat_sample_info_22 <- new_tmt_sample_info[new_tmt_sample_info$Group.b=="5"|new_tmt_sample_info$Group.b=="6"|new_tmt_sample_info$Group.b=="1",]
# F6_127N F5_127N
temp_heat_sample_info_22 <- temp_heat_sample_info_22[!(temp_heat_sample_info_22$MS_ID=="F6_127N"|temp_heat_sample_info_22$MS_ID=="F5_127N"),]

ordered.item12 <- c("5","5","5","5","5","5","5","5","5","5","5","5","5","5","6","6","6","6","6","6","6","6","6","6","6","6","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1")

ordered.item12 <- factor(1:length(ordered.item12),labels = ordered.item12)

temp_heat_sample_info_22$Group.b <- factor(temp_heat_sample_info_22$Group.b,levels = levels(ordered.item12))
temp_heat_sample_info_22 <- temp_heat_sample_info_22[order(temp_heat_sample_info_22$Group.b),]


df4_G45_diff_protein_561 <- df4_G45_diff_protein[,temp_heat_sample_info_22$MS_ID]

new_561_mead_data <- data.frame()
repeat_1 <- unique(temp_heat_sample_info_22$Group.b) ### 2021.12.20
for (i in 1:length(repeat_1)){
  print(i)
  tmp_info <- temp_heat_sample_info_22[temp_heat_sample_info_22$Group.b==repeat_1[i],]
  tmp_ratio <- df4_G45_diff_protein_561[,tmp_info$MS_ID]
  tmp_ratio <- as.data.frame(tmp_ratio)
  tmp_ratio$mean <- apply(tmp_ratio,1,function(x){mean(x)})
  colnames(tmp_ratio) <- c(colnames(tmp_ratio)[1:(dim(tmp_ratio)[2]-1)],as.character(repeat_1[i]))
  if (dim(new_561_mead_data)[1]==0) { new_561_mead_data <- tmp_ratio}
  else {new_561_mead_data <- cbind(new_561_mead_data,tmp_ratio)}
}

new_561_mead_data_2 <- new_561_mead_data[,c("5","6","1")]

new_561_mead_data_2_scale <- t(scale(t(new_561_mead_data_2)))
new_561_mead_data_2_scale_melt <- melt(new_561_mead_data_2_scale)
colnames(new_561_mead_data_2_scale_melt) <- c("protein","Group","Normalized_protein_intensity")
new_561_mead_data_2_scale_melt$Group <- gsub("5","Pre-vaccination",new_561_mead_data_2_scale_melt$Group)
new_561_mead_data_2_scale_melt$Group <- gsub("6","Post-vaccination",new_561_mead_data_2_scale_melt$Group)
new_561_mead_data_2_scale_melt$Group <- gsub("1","Vaccinated-Omicron",new_561_mead_data_2_scale_melt$Group)

new_561_mead_data_2_scale_melt$Group <- factor(new_561_mead_data_2_scale_melt$Group, levels=c("Pre-vaccination","Post-vaccination","Vaccinated-Omicron"))

new_561_line_plot <- ggplot(new_561_mead_data_2_scale_melt,aes(x=Group , y=Normalized_protein_intensity, group=protein)) +
  geom_line(color="gray90",size=0.8) +
  stat_summary(aes(group=1),fun.y=median, geom="line", size=1.2, color="#c51b7d") +
  ylab("Normalized protein intensity")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=8, face = "bold"),
        strip.text = element_text(size = 10, face = "bold"))+
  scale_y_continuous(expand=c(0.01,0.01))+
  scale_x_discrete(expand=c(0.01,0.01))

ggsave("MainFig C.pdf",plot = new_561_line_plot,width = 5,height = 5) 

### boxplot 2022.2.28
# MainFig D
del_sample1 <- c("om1","om12","om1_rep")
new_temp_sample_info <- new_tmt_sample_info
new_temp_sample_info <- new_temp_sample_info[new_temp_sample_info$SampleID %notin% del_sample1,]

# Boxplot 做这几个蛋白(SAA1, SAA2, TGFB1, GPT, CRP, GOT1)在group 1-5的表达情况.
# "P0DJI8_SAA1" "P0DJI9_SAA2" "P01137_TGFB1" "P24298_GPT" "P02741_CRP" "P17174_GOT1"
# df4_plasma_TMT_protein
# boxplot_protein <- c("P0DJI8_SAA1","P0DJI9_SAA2","P01137_TGFB1","P24298_GPT","P02741_CRP","P17174_GOT1")
boxplot_protein <- c("P0DJI8_SAA1","P0DJI9_SAA2","P24298_GPT","P02741_CRP","P17174_GOT1")
df5_rbe_751_slected_protein <- df4_plasma_TMT_protein[boxplot_protein,]

# new_temp_sample_info 图4就用那个boxplot就用删掉om 1 12这两个病人的数据做哈
df5_rbe_751_slected_protein <- df5_rbe_751_slected_protein[,new_temp_sample_info$MS_ID]

df5_rbe_751_slected_protein_boxplot <- melt(df5_rbe_751_slected_protein)
colnames(df5_rbe_751_slected_protein_boxplot) <- c("protein","MS_ID","intensity")
df5_rbe_751_slected_protein_boxplot$Group <- new_tmt_sample_info$Group.b[match(df5_rbe_751_slected_protein_boxplot$MS_ID,new_tmt_sample_info$MS_ID)]

df5_rbe_751_slected_protein_boxplot_2 <- df5_rbe_751_slected_protein_boxplot
df5_rbe_751_slected_protein_boxplot_2$Group <-gsub("3","3/4",df5_rbe_751_slected_protein_boxplot_2$Group)
df5_rbe_751_slected_protein_boxplot_2$Group <-gsub("^4","3/4",df5_rbe_751_slected_protein_boxplot_2$Group)

# PrV PoV OM NO RE
df5_rbe_751_slected_protein_boxplot_2$Group <- gsub("5","PrV",df5_rbe_751_slected_protein_boxplot_2$Group)
df5_rbe_751_slected_protein_boxplot_2$Group <- gsub("6","PoV",df5_rbe_751_slected_protein_boxplot_2$Group)
df5_rbe_751_slected_protein_boxplot_2$Group <- gsub("1","OM",df5_rbe_751_slected_protein_boxplot_2$Group)
df5_rbe_751_slected_protein_boxplot_2$Group <- gsub("2","NO",df5_rbe_751_slected_protein_boxplot_2$Group)
df5_rbe_751_slected_protein_boxplot_2$Group <- gsub("3/4","RE",df5_rbe_751_slected_protein_boxplot_2$Group)


df5_rbe_751_slected_protein_boxplot_2$Group <- factor(df5_rbe_751_slected_protein_boxplot_2$Group, levels=c("PrV","PoV","OM","NO","RE"))

# c("P0DJI8_SAA1","P0DJI9_SAA2","P01137_TGFB1","P24298_GPT","P02741_CRP","P17174_GOT1")
for (protein in boxplot_protein){
  
  temp_df5_rbe_751_slected_protein_boxplot_2 <- df5_rbe_751_slected_protein_boxplot_2[df5_rbe_751_slected_protein_boxplot_2$protein==protein,]
  
  group=levels(factor(temp_df5_rbe_751_slected_protein_boxplot_2$Group))
  temp_df5_rbe_751_slected_protein_boxplot_2$Group=factor(temp_df5_rbe_751_slected_protein_boxplot_2$Group, levels=group)
  comp=combn(group,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  
  p1 <- ggboxplot(temp_df5_rbe_751_slected_protein_boxplot_2, x="Group", y="intensity", color = "Group",  short.panel.labs = FALSE, #### ,palette = "jco"palette = c("#00AFBB", "#E7B800", "#FC4E07") 颜色
                  palette = c("#559FD7", "#143E64", "#A54643","#957CB7","#E9975F"), ###颜色
                  add.params = list(size=0.2),bxp.errorbar = TRUE,bxp.errorbar.width = 0.5,outlier.size=0.02)+
    stat_compare_means(label = "p.format",method ="t.test",comparisons = my_comparisons )+
    xlab("Group")+
    ylab("Protein ratio")
  
  
  ggsave(paste0("MainFig D_",protein,"_boxplot.pdf"),plot =p1 ,device = NULL,width = 5, height = 5)
  
}


pca_SampleID <-  new_tmt_sample_info$Diagnosis[match(colnames(df4_plasma_TMT_protein),new_tmt_sample_info$MS_ID)]
batch_type <-  sapply(colnames(df4_plasma_TMT_protein), function(x) str_split(x,"_")[[1]][1])

# FigS1 A
plot.pca(df4_plasma_TMT_protein,type = batch_type,label2 =NULL,title="FigS1 A",get.color(6),ellise_type = NULL)

df8 <- df4_plasma_TMT_protein
colnames(df8) <- new_tmt_sample_info$SampleID[match(colnames(df8),new_tmt_sample_info$MS_ID)]
pca_SampleID2 <-  new_tmt_sample_info$Diagnosis[match(colnames(df8),new_tmt_sample_info$SampleID)]
batch_type2 <- new_tmt_sample_info$BatchNumber[match(colnames(df8),new_tmt_sample_info$SampleID)]
batch_type2 <- gsub("_","",batch_type2)

# MainFig A
pca_SampleID3 <-  new_tmt_sample_info$Diagnosis._1[match(colnames(df8),new_tmt_sample_info$SampleID)]
plot.pca(df8,type = pca_SampleID3,label2 =NULL,title="MainFig A",get.color(8),ellise_type = NULL)


#########
# df4_plasma_TMT_protein   seqKNN 填充缺失值
# df5_rbe_751   0 填充缺失值
group1 <- new_tmt_sample_info[new_tmt_sample_info$Group.b=="1",]$MS_ID
group3 <- new_tmt_sample_info[new_tmt_sample_info$Group.b=="3",]$MS_ID
group4 <- new_tmt_sample_info[new_tmt_sample_info$Group.b=="4",]$MS_ID

df1_G134 <- df4_plasma_TMT_protein[,c(group1,group3,group4)]
df1_G134 <- as.data.frame(df1_G134)

GO_ID <- read.xlsx("./data/GO_info.xlsx",rowNames = TRUE)
# df1_G134
df_GO <- df1_G134
rownames(df_GO) <- sapply(rownames(df_GO), function(x) str_split(x,";")[[1]][1])
go_id1 <- c()
grep(paste0("\\^","C4B"),rownames(df_GO),value = TRUE)

for (i in 1:dim(GO_ID)[1]) {
  temp_id1 <- grep(paste0("_",rownames(GO_ID)[i],"$"),rownames(df_GO),value = TRUE)
  go_id1 <- c(go_id1,temp_id1)
}

##### 2022.2.28 热图   MainFig B
df2_G134_heat_data_GO <- df1_G134[go_id1,]
rownames(df2_G134_heat_data_GO) <-  gsub("\\^","_",rownames(df2_G134_heat_data_GO))
df2_G134_heat_data_GO_scale <- t(scale(t(df2_G134_heat_data_GO)))
df2_G134_heat_data_GO_scale <- as.matrix(df2_G134_heat_data_GO_scale)

temp_heat_sample_info_3 <- new_tmt_sample_info[new_tmt_sample_info$Group.b=="1"|new_tmt_sample_info$Group.b=="3"|new_tmt_sample_info$Group.b=="4",]
temp_heat_sample_info_3 <- temp_heat_sample_info_3[!(temp_heat_sample_info_3$BatchNumber=="F5_"|temp_heat_sample_info_3$BatchNumber=="F6_"),]

ordered.item112 <- c("om1","om1_rep","om2","om3","om4","om5","om6","om7","om8","om9","om10","om11","om12","om13","om14","om15","om16","om17","flu1","flu2","flu3","flu4","flu5","flu6","flu7","flu8","flu9","flu10","flu11","flu11_rep","flu12","flu13","flu14","flu15","flu16","flu17","flu18","flu19","flu20","ad1","ad2","ad3","ad4","ad5","ad6","ad7","ad8","ad9","ad9_rep")

ordered.item112 <- factor(1:length(ordered.item112),labels = ordered.item112)

temp_heat_sample_info_3$SampleID <- factor(temp_heat_sample_info_3$SampleID,levels = levels(ordered.item112))
temp_heat_sample_info_3 <- temp_heat_sample_info_3[order(temp_heat_sample_info_3$SampleID),]


df2_G134_heat_data_GO_scale <- df2_G134_heat_data_GO_scale[,temp_heat_sample_info_3$MS_ID]

colnames(df2_G134_heat_data_GO_scale) <- new_tmt_sample_info$SampleID[match(colnames(df2_G134_heat_data_GO_scale),new_tmt_sample_info$MS_ID)]

haa = HeatmapAnnotation(Group = as.character(new_tmt_sample_info$Group.b[match(colnames(df2_G134_heat_data_GO_scale),new_tmt_sample_info$SampleID)]),
                        col = list(Group=c("1"=paired_col[6],"2"=paired_col[7],"3"=paired_col[8],"4"=paired_col[9],"5"=paired_col[10],"6"=paired_col[11])))

rowHaa <- rowAnnotation(`Complement System` = GO_ID$Complement.System[match(sapply(rownames(df2_G134_heat_data_GO_scale), function(x) str_split(x,"_")[[1]][2]),rownames(GO_ID))],
                        `Acute.Phase.Response.Signaling` = GO_ID$Acute.Phase.Response.Signaling[match(sapply(rownames(df2_G134_heat_data_GO_scale), function(x) str_split(x,"_")[[1]][2]),rownames(GO_ID))],
                        `NRF2-mediated.Oxidative.Stress.Response` = GO_ID$`NRF2-mediated.Oxidative.Stress.Response`[match(sapply(rownames(df2_G134_heat_data_GO_scale), function(x) str_split(x,"_")[[1]][2]),rownames(GO_ID))],
                        col = list(`Complement System`=c("Complement System"="#6BB03B"),
                                   `Acute.Phase.Response.Signaling`=c("Acute Phase Response Signaling"="#AAC886"),
                                   `NRF2-mediated.Oxidative.Stress.Response`=c("NRF2-mediated Oxidative Stress Response"="#E3EDD8")),
                        na_col = "white")

col_fun=colorRamp2(c(-3:-1,0,1:3), c(brewer.pal(11,"RdYlBu")[11:9],"white",brewer.pal(11,"RdYlBu")[3:1]))
heatmap_11 <- Heatmap(df2_G134_heat_data_GO_scale,
                      name = "Normalized protein ratio",
                      col = col_fun ,
                      na_col = "#CCCCCC",
                      cluster_columns = FALSE,
                      heatmap_height = unit(30,"cm"),
                      heatmap_width  = unit(25,"cm"),
                      top_annotation = haa,
                      right_annotation = rowHaa,
                      show_column_dend = TRUE,
                      show_row_dend = TRUE,
                      show_row_names = TRUE,
                      row_names_gp = gpar(fontsize = 8),
                      show_column_names  = TRUE)

pdf('MainFig B.pdf',width = 15,heigh=15) #保存为pdf文件
draw(heatmap_11)
dev.off()


####### 2022.2.28 pathway enrichment 附件2  
# 气泡图 bubble plot FigS2 A
bubble_plot_data <- read_xls("./data/Group1_34_pathway.xls")
colnames(bubble_plot_data) <- bubble_plot_data[1,]
bubble_plot_data <- bubble_plot_data[2:dim(bubble_plot_data)[1],]

bubble_plot_data_top10 <- bubble_plot_data[1:10,]
bubble_plot_data_top10$`z-score1` <- bubble_plot_data_top10$`z-score`
bubble_plot_data_top10[is.na(bubble_plot_data_top10)] <- c("0.25")
bubble_plot_data_top10$`z-score` <- as.numeric(bubble_plot_data_top10$`z-score`)
bubble_plot_data_top10$`-log(p-value)` <- as.numeric(bubble_plot_data_top10$`-log(p-value)`)

bubble_plot_data_top10 <- bubble_plot_data_top10[order(bubble_plot_data_top10$`-log(p-value)`,decreasing = FALSE),]

pathway <- bubble_plot_data_top10$`Ingenuity Canonical Pathways`
pathway <- factor(1:length(pathway),labels = pathway)
bubble_plot_data_top10$`Ingenuity Canonical Pathways` <- factor(bubble_plot_data_top10$`Ingenuity Canonical Pathways`,levels = levels(pathway))
bubble_plot_data_top10 <- bubble_plot_data_top10[order(bubble_plot_data_top10$`Ingenuity Canonical Pathways`),]

p = ggplot(bubble_plot_data_top10,aes(`-log(p-value)`,`Ingenuity Canonical Pathways`))
p=p + geom_point(aes(size=`-log(p-value)`,fill=`z-score`), colour="black",shape=21)
pbubble = p+ geom_point(aes(size=`-log(p-value)`),shape=21)+scale_size(range = c(5, 12))
pr = pbubble+scale_fill_gradient2(limits = c(-3.2,3),breaks = c(-2,-1,0,1,2),low="#5A8388",high="#A53A34",mid = "white",midpoint = 0.25)

pr = pr+labs(color=expression(z-score),size="-log(p-value)", x="-log(p-value)",y="Pathway name",title="Pathway enrichment")
pp =pr + theme_bw()

ggsave("FigS2 A.pdf",pp, width=3.5, height=2,scale=3)###改变图片比例，大小

# 柱状图 pathway enrichment 附件3 FigS2 B
bar_plot_data_132_1 <- read_xlsx("./data/pathway_group16_56_20220227-2355.xlsx",sheet = "132proteins")

pp_132 <- ggplot(bar_plot_data_132_1, aes(y = reorder(`term description`,-log10(`false discovery rate`))))+ 
  geom_bar(aes(x = -log10(`false discovery rate`),fill= -log10(`false discovery rate`)), stat="identity",colour="black",size=.2,show.legend = TRUE)+
  theme_bw()+
  xlab("−log10(False discovery rate)")+
  ylab("Pathway")+
  theme(
    axis.text.x = element_text(colour="black"),
    axis.text.y = element_text(colour="black"),
    axis.title.x = element_text(colour="black", size=13),
    axis.title.x.top = element_text(colour="black", size=13))+
  scale_fill_continuous(low="#EABED7",high="#CF4D92")

ggsave("FigS2 B.pdf",plot =pp_132 ,device = NULL,width = 10, height = 6)

# FigS2 C
bar_plot_data_333_1 <- read_xlsx("./data/pathway_group16_56_20220227-2355.xlsx",sheet = "333proteins")

pp_333 <- ggplot(bar_plot_data_333_1, aes(y = reorder(`term description`,-log10(`false discovery rate`))))+ 
  geom_bar(aes(x = -log10(`false discovery rate`),fill= -log10(`false discovery rate`)), stat="identity",colour="black",size=.2,show.legend = TRUE)+
  theme_bw()+
  xlab("−log10(False discovery rate)")+
  ylab("Pathway")+
  theme(
    axis.text.x = element_text(colour="black"),
    axis.text.y = element_text(colour="black"),
    axis.title.x = element_text(colour="black", size=13),
    axis.title.x.top = element_text(colour="black", size=13))+
  scale_fill_continuous(low="#F7ECB1",high="#ED9C00")

ggsave("FigS2 C.pdf",plot =pp_333 ,device = NULL,width = 10, height = 6)



