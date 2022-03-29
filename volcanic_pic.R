setwd('D:\\09_纳米水稻\\02_DEG_ss')
data<-read.table('NM_count.txt',header = 1)
library(DESeq2)
library(ggplot2)
countdata<-data[,7:ncol(data)][1:4]
rownames(countdata)<-data[,1]
colnames(countdata)

#表达矩阵数据校正
exprSet <- countdata
boxplot(exprSet,outline=FALSE, notch=T, las=2)
library(limma)
exprSet=data.frame(normalizeBetweenArrays(exprSet))
newdata<-round(exprSet)
boxplot(newdata,outline=FALSE, notch=T, las=2)
# write.table(newdata,file = "D:\\6烟草\\1fpkm_cor\\矫正后的cor/adjust_data.txt",sep = '\t',row.names = T,col.names = T,quote = F)

#判断数据是否需要转换
# 1,2     BCD
# 3,4     CON
# 5,6     YCD
#                                 DEseq
colnames(newdata)
group_list=factor(c(rep('BCD',2),rep('CK',2)))
# group_list <- factor(group_list,levels = c("Normal","Tumor"),ordered = F)
exprSet<-newdata[,c(c(1,2),c(3,4))]
# exprSet <- log2(exprSet+1)
colData<-data.frame(row.names = colnames(exprSet),
                    group_list=group_list)
colData
colnames(exprSet)
#下面这个为什么叫dds
dds<-DESeqDataSetFromMatrix(countData = exprSet, 
                            colData = colData,  
                            design = ~group_list)     
dds  
# dds <- dds[rowSums(counts(dds))> 10,]    
# dds$group_list <- relevel(dds$group_list, ref = "CON")
dds$group_list <- relevel(dds$group_list, ref = "CK")
dds<-DESeq(dds) 
resultsNames(dds)
res<-results(dds,name = "group_list_BCD_vs_CK") 
#res <- results(dds, contrast=c("group_list","D4","D12")) #返回差异分析结果
res #一开始是一些描述性信息，跟矩阵是不一样的
resOrdered <- res[order(res$padj),] 
head(resOrdered) 
DEG=as.data.frame(resOrdered) 
DEG=na.omit(DEG)
write.table(DEG,file = './group_list_BCD_vs_CK_allDEG.txt',sep = '\t',quote = FALSE)
#..............................差异分析到这里就可以了，差异基因也有了
# res<-read.table("D:/ARyuyan/R_data/L12VSD4.Rdata")
# res<-read.csv("D:/ARyuyan/R_data/group_list_D4_vs_D12/DEG")
# row.names(res)<-res[,1]
# res<-res[,-1]
# DEG<-read.table('D:\\8_wheat_MH\\6_rna-seq/DEG.txt')
diff_gene <- subset(DEG,padj<0.05 )#& abs(log2FoldChange) >=1
# diff_gene <- subset(DEG,padj<0.05) 
diffdata <- as.data.frame(diff_gene)
a<-row.names(diffdata)
length(a)
# write.table(a,file = "D:\\6烟草\\2count_DEG/DEG.csv",row.names = F,col.names = F,quote = F)

gene_up <- row.names(diff_gene[diff_gene$log2FoldChange > 0, ])
length(gene_up)   
up<-data.frame(gene_up)
write.table(up,file = "./group_list_YCD_vs_BCD-up1.txt",row.names = F,col.names = F,quote = F)
gene_down <- row.names(diff_gene[diff_gene$log2FoldChange < 1, ])
length(gene_down) 
down<-data.frame(gene_down)
write.table(down,file = "./group_list_YCD_vs_BCD-down1.txt",row.names = F,col.names = F,quote = F)



#########画图


#logFC_cutoff <- with(DEG,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)) ) 
logFC_cutoff=0 


DEG$change = as.factor(ifelse(DEG$padj < 0.05 & abs(DEG$log2FoldChange) >=logFC_cutoff, 
                              ifelse(DEG$log2FoldChange >=logFC_cutoff ,'Up-regulated','Down-regulated'),'No significant') 
) 
this_tile <- paste0('Cutoff for log2FoldChange is ',round(logFC_cutoff,3),
                    '\nCutoff for Padj is 0.05') 
# '\nThe number of up gene is ',nrow(DEG[DEG$change =='Up-regulated',]) , 
# '\nThe number of down gene is ',nrow(DEG[DEG$change =='Down-regulated',]) 
# ) 

library(ggplot2)
ggplot(data=DEG,  
       aes(x=log2FoldChange, y=-log10(padj),  
           color=change)) + 
  geom_point(alpha=0.6, size=2) + 
  theme_set(theme_set(theme_bw(base_size=20)))+ 
  xlab("log2(Fold change)") + ylab("-log10(Padj)") + 
  # ggtitle( this_tile ) +
  theme(plot.title = element_text(size=15,hjust = 0.5))+ 
  scale_colour_manual(values = c('#20B2AA','gray','#FF9900')) +
  # geom_hline(aes(yintercept=-0.5),colour="#808080",linetype="dashed")+
  # geom_vline(aes(xintercept=1),colour="#808080",linetype="dashed")+
  # geom_vline(aes(xintercept=1),colour="#808080",linetype="dashed")+
  annotate('text', label = (nrow(DEG[DEG$change =='Up-regulated',])), 5,8,size=5.5,fontface  = 'bold',colour = "black") +
  annotate('text', label = (nrow(DEG[DEG$change =='Down-regulated',])), -5, 8,size=5.5,fontface  = 'bold',colour = "black")+
  theme_classic()+
  xlim(-9,9)+ylim(0,14)+
  geom_hline(aes(yintercept=-log10(0.05)),colour='#030303',linetype='dashed')+
  guides(colour = guide_legend(title = 'BCD vs CK',override.aes = list(size=3, stroke=1.5)))+
  theme(strip.text.x = element_text( colour = 'black',size = 13),#size=8,angle=75
        strip.text.y = element_text( colour = 'black',size = 13),
        strip.background = element_rect(colour="white", fill="grey"),
        axis.title = element_text(face = "bold",
                                  size = "16",color = "black"),
        axis.text.x = element_text(face = "bold",color = "black",
                                   size = 12,  hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(face = "bold",size = 12,color = "black"),
        legend.text = element_text(size = 12,face = "bold"),
        legend.title = element_text(size = 12,
                                    face = "bold"),
        legend.position ='right',
        # legend.key.size = 
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(size=1, colour = "black"),     #####################     坐标轴边框加粗（只能加粗左下轴）
        axis.ticks = element_line(size = 1),
        # panel.border = element_rect(colour = "black", fill=NA, size=1),
        # axis.line = element_line(size=2, colour = "black"),
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)) 

# ggsave(g,filename = 'D:/ARyuyan/R_data/Anewstandspecific/count/group_list_L12_vs_D4/volcano.png') 






a<-row.names(DEG)
write.table(a,file ="D:/ARyuyan/R_data/strand-specific/DEGgroup_list_D12rev_vs_D12for.csv",row.names =F,col.names = T,quote = F,sep = "\t" )
#画个热图
DEG=na.omit(DEG)  #去除一些NA的值

nrDEG=DEG       #要改一下=后后边的数据
## heatmap 热图
library(pheatmap) 
choose_gene=head(rownames(nrDEG),50) ## 50 maybe better 
choose_matrix=exprSet[choose_gene,] 
choose_matrix=t(scale(t(choose_matrix))) 
pheatmap(choose_matrix,filename = 'D:/ARyuyan/R_data/group_list_L4_vs_L12/2group_list_L4_vs_L12.png') 

## volcano plot 
colnames(nrDEG) 
#plot(nrDEG$logFC,-log10(nrDEG$pvalue)) 


DEG=nrDEG 


#。。。............................edgeR

library(edgeR) 
d <- DGEList(counts=exprSet,group=factor(group_list)) 
length(d)  #看d的元素，这里d包括了表达矩阵和分组信息
d$samples$lib.size <- colSums(d$counts) #算了一下reads的大小
d <- calcNormFactors(d) #计算了Normalization的因子，根据某个比例进行缩放他
d$samples #normalization的大小接近于1说明差不多了

#这个时候需要一个d矩阵
dge=d #因为后续代码都为dge所以转化一下

design <- model.matrix(~0+factor(group_list)) 
rownames(design)<-colnames(dge) 
colnames(design)<-levels(factor(group_list)) 

dge <- estimateGLMCommonDisp(dge,design) #以下三个dge要看统计学原理才能看到懂，
#如果只要做差异分析就运行一下就连可以了，之所以要运行三次，因为有时
#候不只是两两比较，还会有更多矩阵比较，这个时候需要你对你的数据充分理解
dge <- estimateGLMTrendedDisp(dge, design) 
dge <- estimateGLMTagwiseDisp(dge, design) 

fit <- glmFit(dge, design) 

lrt <- glmLRT(fit,  contrast=c(1,0)) #contrast=c(-1,1,0))修改要看design
nrDEG=topTags(lrt, n=nrow(exprSet)) 
nrDEG=as.data.frame(nrDEG) #差异分析结果也有了，有两万多个
head(nrDEG) 

#。。。。。。。。。。。。。。到这里差异分析就结束了
write.csv(nrDEG,"DEG_treat_12_edgeR.csv",quote = F) 

#lrt <- glmLRT(fit, contrast=c(-1,0,1) ) 
# nrDEG=topTags(lrt, n=nrow(exprSet)) 
#  nrDEG=as.data.frame(nrDEG) 
# head(nrDEG) 
# write.csv(nrDEG,"DEG_treat_2_edgeR.csv",quote = F) 
# summary(decideTests(lrt)) 
#  plotMD(lrt) 
#  abline(h=c(-1, 1), col="blue") 
# } 

#。。。。。。。。。。。。。。。limma
BiocManager::install('limma')
suppressMessages(library(limma)) 
design <- model.matrix(~0+factor(group_list)) 
colnames(design)=levels(factor(group_list)) 
rownames(design)=colnames(exprSet) 


dge <- DGEList(counts=exprSet) 
dge <- calcNormFactors(dge) 
logCPM <- cpm(dge, log=TRUE, prior.count=3) #这个logcpm值没什么用，因为在后面没用到


v <- voom(dge,design,plot=TRUE, normalize="quantile") #voom直接根据dge，design来算
fit <- lmFit(v, design) #这个V就可以根据design来作出差异分析


group_list
#下面需要修改一下，要知道哪个和哪个比较
cont.matrix=makeContrasts(contrasts=c('D4-L4'),levels = design) 

fit2=contrasts.fit(fit,cont.matrix) 
fit2=eBayes(fit2) 

tempOutput = topTable(fit2, coef='D4-L4', n=Inf) 
DEG_limma_voom = na.omit(tempOutput) 

#画个热图
DEG=na.omit(DEG_limma_voom)  #去除一些NA的值
nrDEG=DEG       #要改一下=后后边的数据
## heatmap 热图
library(pheatmap) 
choose_gene=head(rownames(nrDEG),50) ## 50 maybe better 
choose_matrix=exprSet[choose_gene,] 
choose_matrix=t(scale(t(choose_matrix))) 
pheatmap(choose_matrix,filename = 'DEG_limma_voomheatmap.png') 




