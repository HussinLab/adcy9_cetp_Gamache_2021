#!/usr/bin/env Rscript

library('ggplot2')
print('Figure 5d')

######################################
#####		 	  				######
#####		 READ  				######
#####		 	  				######
######################################

pdf('CETP.rs1967309-rs158477.Hg38.GTEx.SkinExposed.25PEER.pdf',width=7, height=7)


# This file contains residual of CETP expression made on CETP expression, corrected for all covariates (no SNPs)
allData=read.table('GTEx.Skin.Male.25PEERs.txt',header=TRUE, row.names=NULL, sep=' ')

allData$rs1967309=as.factor(allData$rs1967309)
print(ggplot(allData,aes(y=Resid,fill=rs1967309,x=Geno))+geom_boxplot(outlier.shape=NA,color='black')+geom_jitter()+theme(panel.background=element_rect(fill='white', color='black'), legend.position='none')+scale_fill_manual(values=c('#56bb5d','#317d36','#0b3c0e'))+ylim(-4,3.5))



dev.off()
