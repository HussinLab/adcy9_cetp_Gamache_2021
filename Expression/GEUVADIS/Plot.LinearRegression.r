#!/usr/bin/env Rscript

library('ggplot2')
print('Figure 4b')

######################################
#####		 	  				######
#####		 READ  				######
#####		 	  				######
######################################

pdf('CETP.European.rs1967309-rs158477.Hg38.1000G.PEER.pdf',width=7, height=7)


# This file contains residual of CETP expression made on CETP expression, corrected for all covariates (no SNPs)
allData=read.table('ResidualExpressionCETP.txt',header=TRUE, row.names=NULL, sep='\t')

# Order genotype to be in good order in graph
allData$Geno=factor(allData$Geno,c('AA GG','AA GA','AA AA','AG GG','AG GA','AG AA','GG GG','GG GA','GG AA'))
allData$rs158477_Geno=factor(allData$rs158477_Geno,c('GG','GA','AA'))


allData$rs1967309=as.factor(allData$rs1967309)
print(ggplot(allData,aes(y=resid,fill=rs1967309,x=Geno))+geom_boxplot(outlier.shape=NA,color='black')+geom_jitter()+theme(panel.background=element_rect(fill='white', color='black'), legend.position='none')+scale_fill_manual(values=c('#56bb5d','#317d36','#0b3c0e'))+ylim(-4,3.5))
        # +theme(legend.position = "none",
        # panel.grid = element_blank(),
        # axis.title = element_blank(),
        # axis.text = element_blank(),
        # axis.ticks.y = element_blank()))



dev.off()
