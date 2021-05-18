library(ggplot2)
library('reshape2')

print('Figure 5a')

png('UKphewas.Interaction.png', height=400,width=1200)


data=read.table('PheWas.Interaction.txt',header=TRUE)
data=data[order(-data$type3ss_pvalue_d_CETP_d_ADCY9),]

data$Outcome_name=factor(data$Outcome_name, levels=unique(data$Outcome_name))

ggplot(data, aes(x=Outcome_name, y=type3ss_pvalue_d_CETP_d_ADCY9, color=sex2))+geom_point(size=3)+theme(axis.text.x=element_text(angle=90, vjust=0.5,hjust=1, color='black'), legend.position='none', panel.background=element_rect(fill=NA), panel.border=element_rect(colour='black', fill=NA), axis.text=element_text(size=20))+geom_hline(yintercept=c(-log10(0.10),-log10(0.05)),col=c('orange','red'), linetype='dashed')+scale_color_manual(values=c('#E69F00','#009E63','#CC79A7'))



