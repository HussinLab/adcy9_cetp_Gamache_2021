#ADCY9 Project 
#only KO samples
#RSEM rf-stranded

options(width=10000)

library(limma)
library(edgeR)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(ggfortify)
library(ggrepel)
library(qvalue)
library(sva)
library(pheatmap)
library(RColorBrewer)
library(dendsort)

get_contrast=function(fit,design,contrast,name,table){
    vec=rep(0,length(colnames(design)))
    vec[abs(contrast)]=contrast/abs(contrast)
    fit2 <- contrasts.fit(fit, vec)
    fit2 <- eBayes(fit2)
    
    tab=data.frame(beta=fit2$coefficients,p=fit2$p.value)
    tab$fdr=qvalue(tab$p)$qvalues
    
    colnames(tab)=paste0(colnames(tab),"_",name)
    
    print(paste(name,length(which(tab[,3]<0.01)),length(which(tab[,3]<0.05)),length(which(tab[,3]<0.1)),length(which(tab[,3]<0.2))))
    
    if(missing(table)){
        return(tab)
    }else{
        tab=cbind(table,tab)
        return(tab)
    }
}



cols=read.table("df.txt",header=TRUE,sep="\t")
row.names(cols)=cols$SampleID
reads=read.table("ADCY9_Tardif.KO.GRCh38.E87.no205.genes.txt",header=TRUE,row.names=1)



reads= reads[,colnames(reads) %in% cols$SampleID]
reads= reads[,order(colnames(reads))]


gene_infos=read.table("Homo_sapiens.GRCh38.87.genes.tsv",header=T)
coding_ids=read.table("Homo_sapiens.GRCh38.87.protein_coding.genes.ids.txt",header=F)
coding_ids_hugo = merge(coding_ids,gene_infos,by.x="V1",by.y="gene_id")

countsData=reads[rownames(reads) %in% coding_ids$V1,]

system("mkdir KO_filter10mean")



#OVERALL PCA



cols=cols[order(cols$SampleID),]

counts_tmp= countsData[,colnames(countsData) %in% cols$SampleID]
counts_tmp= counts_tmp[,order(colnames(counts_tmp))]



#QCs

dge <- DGEList(
    samples = cols,
    counts = counts_tmp
    )



# Use the conditions as groups
dge$samples$group <- dge$samples$Condition

# Convert TPM values to log cpm
lcpm <- cpm(dge, log=TRUE)

# Median of expression
rowMedians <- function(x, na.rm = FALSE) {
    apply(x, 1, median, na.rm = na.rm)
}

melt_lcpm <- melt(lcpm)

median_melt_lcpm <- melt_lcpm %>%
  group_by(Var1) %>%
  summarize(value = median(value))

## Compare each distribution against the median values accross the whole dataset.
## Makes spotting bad samples easier
system("mkdir -p Output/Figures")

pdf('Output/Figures/LCPM_distribution.density.pdf', width = 14)
p <- ggplot(data = melt_lcpm, aes(x = value)) + facet_wrap(~ Var2) + geom_density(color = 'deepskyblue3') + geom_density(data = median_melt_lcpm, aes(x = value))+theme()
print(p)
dev.off()



## Thresholding: filter out lowly expressed genes that could throw off the analyses
# Keep genes for which the mean lcpm value is above 0
average_lcpm <- aveLogCPM(dge)

pdf('Output/Figures/Average_LCPM_distribution.density.pdf')
p <- ggplot(data = data.frame(val = average_lcpm), aes(x = val)) + geom_density()
print(p)
dev.off()

# Table of genes with an average log2 count per million above 0
summary(average_lcpm > 0)
#   Mode   FALSE    TRUE
#logical    8066   11895

completeGeneList=row.names(dge)
#test another threshold
passingLCPM=completeGeneList[average_lcpm > 0]

toExtract=intersect(coding_ids$V1,passingLCPM)

dge <- dge[toExtract , , keep.lib.sizes=FALSE]

# Show the effect of normalization
raw_dge <- dge
norm_dge <- calcNormFactors(dge, method = "TMM")

raw_lcpm <- cpm(raw_dge, log = TRUE)
norm_lcpm <- cpm(norm_dge, log = TRUE)

df_new <- rbind(
    {
        k <- melt(raw_lcpm)
        k$group = factor('raw')
        k
    },
    {
        k <- melt(norm_lcpm)
        k$group = factor('normalized')
        k
    }
    )

# Boxplots showing the distribution of expressions for each individual, before and after normalization
pdf('Output/Figures/LCPM_normalization.boxplot_LCPM_2.pdf', width = 24)
p <- ggplot(data = df_new, aes(x = Var2, y = value, color = group)) + geom_boxplot() + facet_wrap(~group, ncol=1) + guides(color = FALSE) + theme(axis.text.x=element_text(angle=90,hjust=1))
print(p)
dev.off()



sv_result <- svd(norm_dge)

df <- data.frame(
    x = seq(1:length(sv_result$d)),
    y = sv_result$d^2/sum(sv_result$d^2)
    )

pdf('Output/Figures/Explained_variance_LCPM_2.pdf')

p <- ggplot(data = df, aes(x = x, y = y)) + 
    geom_point() + 
    geom_line(linetype = 'dashed') + 
    coord_cartesian(xlim = c(1, 10)) +
    labs(title = 'Standard Value Decomposition', x = 'Principal component', y = 'Proportion of explained variance')
print(p)

dev.off()



#KO analysis

samples_ko = rownames(subset(cols,Type=="KO" ))

sampleTable=cols[samples_ko,]

sampleTable.bkp=sampleTable
sampleTable$SampleID=factor(sampleTable$SampleID,levels(unique(as.factor(as.character(sampleTable$SampleID)))))
sampleTable$Type=factor(sampleTable$Type,levels(unique(as.factor(as.character(sampleTable$Type)))))
sampleTable$Condition=factor(as.factor(sampleTable$Condition),levels(unique(as.factor(as.character(sampleTable$Condition)))))




counts_tmp= countsData[,colnames(countsData) %in% sampleTable$SampleID]
counts_tmp= counts_tmp[,order(colnames(counts_tmp))]


dge <- DGEList(
    samples = sampleTable,
    counts = counts_tmp
    )



# Use the conditions as groups
dge$samples$group <- dge$samples$Condition

# Convert TPM values to log cpm
lcpm <- cpm(dge, log=TRUE)

# Median of expression
rowMedians <- function(x, na.rm = FALSE) {
    apply(x, 1, median, na.rm = na.rm)
}

melt_lcpm <- melt(lcpm)

median_melt_lcpm <- melt_lcpm %>%
  group_by(Var1) %>%
  summarize(value = median(value))


## Thresholding: filter out lowly expressed genes that could throw off the analyses
# Keep genes for which the mean lcpm value is above 0
average_lcpm <- aveLogCPM(dge)

# Table of genes with an average log2 count per million above 1
summary(average_lcpm > 0)
#   Mode   FALSE    TRUE
#logical    8168   11793


#TRY ANOTHER ONE....
passingMean10 = data.frame(Control=rowMeans(counts_tmp[,subset(sampleTable,Condition=="Control")$SampleID])>10 , Alt1=rowMeans(counts_tmp[,subset(sampleTable,Condition=="Alt1")$SampleID])>10)

summary(rowSums(passingMean10)>0)
   # Mode   FALSE    TRUE
# logical    6795   13166

completeGeneList=row.names(dge)
#test another threshold
# passingLCPM=completeGeneList[average_lcpm > 0]

passingMean10_1group=completeGeneList[rowSums(passingMean10)>0]

# toExtract_ko=intersect(coding_ids$ENSEMBL_ID,passingLCPM)
toExtract_ko=intersect(coding_ids$V1,passingMean10_1group)





dge <- dge[toExtract_ko , , keep.lib.sizes=FALSE]

# Use the conditions as groups
dge$samples$group <- dge$samples$Condition


raw_dge <- dge
norm_dge <- calcNormFactors(dge, method = "TMM")

raw_lcpm <- cpm(raw_dge, log = TRUE)
norm_lcpm <- cpm(norm_dge, log = TRUE)
norm_cpm <- cpm(norm_dge, log = FALSE)



#PCA
#mat_colors_tmp <- list(group = brewer.pal(3, "Set1"))
#names(mat_colors_tmp$group) = c("Control","Alt1","Alt2")
mat_colors_tmp <- list(group = c("red","blue"))
names(mat_colors_tmp$group) = c("Control","Alt1")

sampleTable$Condition = factor(sampleTable$Condition,levels=c("Control","Alt1"))
sampleTable$Pairing=c(1,2,3,4,5,1,2,3,4,5)
sampleTable$Pairing = factor(sampleTable$Pairing,levels=c("1","2","3","4","5"))

pca=prcomp(t(norm_lcpm))
pdf("KO_filter10mean/PCA.uncorrected.normalizedReads.allSamples.pdf")
autoplot(pca,data=sampleTable,colour='Condition',shape='Pairing')+scale_color_manual(values=mat_colors_tmp$group[as.character(sampleTable$Condition)])
dev.off()


#Correct for pairing
new.x= apply(as.matrix(norm_lcpm), 1, FUN = function(x){return(resid(lm(as.numeric(x) ~ +as.factor(sampleTable[,'Pairing']), na.action=na.exclude)))})
row.names(new.x)=colnames(norm_lcpm)


pca=prcomp(new.x)
pdf("KO_filter10mean/PCA.correctedPairing.normalizedReads.allSamples.pdf")
autoplot(pca,data=sampleTable,colour='Condition')+scale_color_manual(values=mat_colors_tmp$group[as.character(sampleTable$Condition)])
dev.off()





y <- dge
y <- calcNormFactors(y)

mod <- model.matrix(~0+Condition+Pairing, data = sampleTable)

design=mod


v <- voom(y,design,plot=FALSE)
fit <-lmFit(v,design)
fit <- eBayes(fit)

voomed_reads = as.data.frame(v$E)

#> data.frame(colnames(design))
# 1 ConditionControl
# 2    ConditionAlt1
# 3         Pairing2
# 4         Pairing3
# 5         Pairing4
# 6         Pairing5

contrast_Alt1_Control_ko=c(2,-1)



results=get_contrast(fit,design,contrast_Alt1_Control_ko,"Alt1_Control_ko")

results=results[,order(colnames(results))]

# Alt1_Control_ko 8400 10457 11593 13166



#Pairing CORRECT FDR 1% + log2fc > 1
dim(subset(results,fdr_Alt1_Control_ko<0.01 & abs(beta_Alt1_Control_ko) > 1 ))
# [1] 1245    9

#Pairing CORRECT FDR 5% + log2fc > 1
dim(subset(results,fdr_Alt1_Control_ko<0.05 & abs(beta_Alt1_Control_ko) > 1 ))
# [1] 1268   9


#write results by FDR thresholds w/o beta
write.table(row.names(subset(results,fdr_Alt1_Control_ko<0.05)),file="KO_filter10mean/Results_Alt1_Control_ko.paired.fdr5perc.txt",row.names=F,quote=F,col.names=F)
write.table(row.names(subset(results,fdr_Alt1_Control_ko<0.01 & abs(beta_Alt1_Control_ko)>1 )),file="KO_filter10mean/Results_Alt1_Control_ko.paired.fdr1perc_absbeta1.txt",row.names=F,quote=F,col.names=F)

#up regulated
write.table(row.names(subset(results,fdr_Alt1_Control_ko<0.05 & beta_Alt1_Control_ko > 0)),file="KO_filter10mean/Results_Alt1_Control_ko.paired.fdr5perc.upregulated.txt",row.names=F,quote=F,col.names=F)

#down regulated
write.table(row.names(subset(results,fdr_Alt1_Control_ko<0.05 & beta_Alt1_Control_ko < 0)),file="KO_filter10mean/Results_Alt1_Control_ko.paired.fdr5perc.downregulated.txt",row.names=F,quote=F,col.names=F)

write.table(results,file="KO_filter10mean/Results_KO.paired.txt",quote=F)

write.table(row.names(results),file="KO_filter10mean/background_KO.txt",row.names=F,quote=F,col.names=F)

#results_lungs_sva9_only2CdtInModel_valsartan_PBI = results

results_ko = results



#DRAW VOLCANO PLOTS
results_ko_withHUGO=merge(results_ko,coding_ids_hugo,by.x="row.names",by.y="V1",all.x=TRUE)
rownames(results_ko_withHUGO)=results_ko_withHUGO$Row.names
results_ko_withHUGO=results_ko_withHUGO[,-1]



volcano=results_ko_withHUGO

pdf("KO_filter10mean/volcano.ko_si-1039_vs_Controls.FDR1perc.l2fc1.notext.pdf")
volcano$threshold <- 
  factor(ifelse(volcano$beta_Alt1_Control_ko>1 & volcano$fdr_Alt1_Control_ko< 0.01, "Upregulated",
                ifelse(volcano$beta_Alt1_Control_ko< -1 & volcano$fdr_Alt1_Control_ko< 0.01,
                       "Downregulated",
                        "Not Significant")
                )
  )

thresholdColors=c("red","darkgreen","blue")
names(thresholdColors)=levels(volcano$threshold)

ggplot(data=volcano, aes(x=beta_Alt1_Control_ko, y=-1*log10(p_Alt1_Control_ko)) ) +
  geom_point(alpha=0.4, size=1,aes(colour=thresholdColors[volcano$threshold])) + ylim(c(0,17))+ggtitle("KO : FDR 1% and abs(log2FC)>1 si-1039 vs SCR") + scale_color_manual(name='Threshold', labels = c('Upregulated','Not Significant','Downregulated'),values=c('blue','darkgreen','red'))



dev.off()




pdf("KO_filter10mean/volcano.ko_si-1039_vs_Controls.FDR1perc.l2fc1.pdf")
volcano$threshold <- 
  factor(ifelse(volcano$beta_Alt1_Control_ko>1 & volcano$fdr_Alt1_Control_ko< 0.01, "Upregulated",
                ifelse(volcano$beta_Alt1_Control_ko< -1 & volcano$fdr_Alt1_Control_ko< 0.01,
                       "Downregulated",
                        "Not Significant")
                )
  )

thresholdColors=c("red","darkgreen","blue")
names(thresholdColors)=levels(volcano$threshold)

ggplot(data=volcano, aes(x=beta_Alt1_Control_ko, y=-1*log10(p_Alt1_Control_ko)) ) +
  geom_point(alpha=0.4, size=1,aes(colour=thresholdColors[volcano$threshold])) + geom_text_repel(
        data=subset(volcano,abs(beta_Alt1_Control_ko)>1 & fdr_Alt1_Control_ko<0.0000000001), aes(label = gene_name),
        box.padding = 0.35, point.padding = 0.5,size=3,
        segment.color = 'grey50')+ ggtitle("KO : FDR 1% and abs(log2FC)>1 si-1039 vs SCR") + scale_color_manual(name='Threshold', labels = c('Upregulated','Not Significant','Downregulated'),values=c('blue','darkgreen','red'))



dev.off()





pdf("KO_filter10mean/volcano.ko_si-1039_vs_Controls.FDR1perc.l2fc1.ADCY9.CETP.pdf")
volcano$threshold <- 
  factor(ifelse(volcano$beta_Alt1_Control_ko>1 & volcano$fdr_Alt1_Control_ko< 0.01, "Upregulated",
                ifelse(volcano$beta_Alt1_Control_ko< -1 & volcano$fdr_Alt1_Control_ko< 0.01,
                       "Downregulated",
                        "Not Significant")
                )
  )

thresholdColors=c("red","darkgreen","blue")
names(thresholdColors)=levels(volcano$threshold)

ggplot(data=volcano, aes(x=beta_Alt1_Control_ko, y=-1*log10(p_Alt1_Control_ko)) ) +
  geom_point(alpha=0.4, size=1,aes(colour=thresholdColors[volcano$threshold])) + geom_text_repel(
        data=subset(volcano,gene_name %in% c("ADCY9","CETP")), aes(label = gene_name),
        box.padding = 0.35, point.padding = 0.5,size=3,
        segment.color = 'grey50')+ ggtitle("KO : FDR 1% and abs(log2FC)>1 si-1039 vs SCR") + scale_color_manual(name='Threshold', labels = c('Upregulated','Not Significant','Downregulated'),values=c('blue','darkgreen','red'))


dev.off()

