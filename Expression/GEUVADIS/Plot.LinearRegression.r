#!/usr/bin/env Rscript

library('ggplot2')
print('Figure 4b')

######################################
#####		 	  				######
#####		 READ  				######
#####		 	  				######
######################################

pdf('CETP.European.rs1967309-rs158477.Hg38.1000G.PEER.pdf',width=5,height=5)
# This file contains residual of CETP expression made on CETP expression, corrected for all covariates (no SNPs)
allData=read.table('ResidualExpressionCETP.txt',header=TRUE, row.names=NULL, sep='\t')


allData$Geno=factor(allData$Geno,c('AA GG','AA GA','AA AA','AG GG','AG GA','AG AA','GG GG','GG GA','GG AA'))
allData$rs158477_Geno=factor(allData$rs158477_Geno,c('GG','GA','AA'))

print(ggplot(allData,aes(y=resid,fill=rs158477_Geno,x=rs158477_Geno))+geom_boxplot(outlier.shape=NA)+geom_jitter()+theme(panel.background=element_rect(fill='white', color='black'), legend.position='none')+scale_fill_manual(values=c('grey30','grey60','grey90')))
allData$rs1967309=as.factor(allData$rs1967309)
print(ggplot(allData,aes(y=resid,fill=rs1967309,x=Geno))+geom_boxplot(outlier.shape=NA)+geom_jitter()+theme(panel.background=element_rect(fill='white', color='black'), legend.position='none')+scale_fill_manual(values=c('#56bb5d','#317d36','#0b3c0e')))



dev.off()
