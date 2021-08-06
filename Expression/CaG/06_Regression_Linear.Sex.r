#!/usr/bin/env Rscript
library('ggplot2')
print('Figure 5c')



###################################################
####                                          #####
####          Variable				          #####
####                                          #####
###################################################



pdf='PEER.number.CaG.Sex.pdf'
pdf(pdf,width=7, height=7)

allValueTissu=readRDS('ENSG00000087237.rds')


print(head(allValueTissu))

print(ggplot(allValueTissu,aes(y=-log10(Pvalue),col=Sex,x=NbPEER, shape=as.character(DirectionBeta)))+geom_line(size=2)+geom_point(aes(fill=Coding), size=5,stroke=2)+scale_color_manual(values=c("#6e4bf8","#b31700"))+scale_fill_manual(values=c("green","orange"))+scale_shape_manual(values=c(unique((allValueTissu$DirectionBeta))))+theme(panel.background=element_rect(fill='white', color='black'), legend.position='none')+geom_hline(lty=2,col=c('orange','red','pink'), lwd=1.5,yintercept=-log10(c(0.1,0.05,0.01))))

# "#6e4bf8" = Male
# "#b31700" = Female
# "green" = additive coding
	# Triangle = positive beta
	# Inversed triangle = negative beta
# "orange" = genotypic coding



	