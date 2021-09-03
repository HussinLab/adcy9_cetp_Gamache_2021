library(ggplot2)
library('reshape2')
library(RColorBrewer)


print('Figure 6-figure supplement 1')

pdf('UKphewas.v1.pdf', height=7,width=13)


data=readRDS('PheWas.NoInteraction.rds')

ggplot(data, aes(x=outcome_label, y=type3ss_pvalue, color=sex2, shape=variant_id))+geom_point(size=3)+theme(axis.text.x=element_text(angle=90, vjust=0.5,hjust=1, color='black'), panel.background=element_rect(fill=NA), panel.border=element_rect(colour='black', fill=NA), axis.text=element_text(size=20))+geom_hline(yintercept=c(-log10(0.05)),col=c('black'), linetype='dashed')+scale_color_manual(values=c('#E69F00','#009E63','#CC79A7'))
