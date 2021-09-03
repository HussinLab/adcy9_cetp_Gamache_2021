#!/usr/bin/env Rscript
library("gplots")
library("RColorBrewer")

print('Figure 6-figure supplement 1b')

pheno='MHHRTATT' 

snp1='rs1967309'
snp2='rs158477'


pdf=paste(c('Combination.',pheno,'.',snp1,'_',snp2,'.pdf'),collapse="")

pdf(pdf,width=5, height=5)

mypalette<-colorRampPalette(c("white", 'cyan', "purple"))(n = 100)
nb_combi=readRDS('Heatmap.x.rds')
nb_ratio=readRDS('Heatmap.cellnote.rds')

heatmap.2(x=nb_combi, Rowv=FALSE, col=mypalette,  Colv=FALSE, dendrogram="none", cellnote=nb_ratio, notecex=1.5, notecol="black", colsep=c(0,1,2,3), rowsep=c(0,1,2,3),sepcol="black",
	sepwidth=c(0.01,0.01), trace="none",margin=c(6,5), symm=F,symkey=F,symbreaks=F, scale="none", ,density.info='none',
	xlab=snp1, ylab=snp2,
	key=FALSE, lwid=c(0.1,0.6), lhei=c(0.1,0.5), cex.lab=1.5,cexRow=1.5,cexCol=1.5)
 



