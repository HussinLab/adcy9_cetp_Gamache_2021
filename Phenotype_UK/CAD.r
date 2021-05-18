library(ggplot2)
library('reshape2')

print('Figure 5b')

# Cardiovascular

data=read.table('CAD.txt', sep='\t', header=TRUE)

# Factorize data
data$rs1967309_fac=factor(data$rs1967309_fac, levels=c('G/G','A/G','A/A'))
data[data$male==1,'male']='Men'
data[data$male==0,'male']='Women'
data$male=factor(data$male, levels=c('Men','Women'))


data$endpoint=factor(data$endpoint, levels=c('cad_prevalent_or_incident','mi_prevalent_or_incident'))

col_7309=c('#56bb5d','#317d36','#0b3c0e','green')
names(col_7309)=c('G/G','A/G','A/A')
colScale=scale_colour_manual(name='rs1967309_fac', values=col_7309)
cad=data







# Biomarkers

data=read.table('BM.txt', sep='\t', header=TRUE)


data$rs1967309_fac=factor(data$rs1967309_fac, levels=c('G/G','A/G','A/A'))
data[data$male==1,'male']='Men'
data[data$male==0,'male']='Women'
data$male=factor(data$male, levels=c('Men','Women'))
data$endpoint=factor(data$endpoint, levels=c('lpa_std','crp_std','hdl_std','apoa_std','ldl_std','apob_std'))

col_7309=c('#56bb5d','#317d36','#0b3c0e','green')
names(col_7309)=c('G/G','A/G','A/A')
colScale=scale_colour_manual(name='rs1967309_fac', values=col_7309)

bm=data


# Plot results

png('UK.CAD.png', height=150,width=1200)
xlim=range(c(exp(cad$AME),exp(cad$AME-cad$SE*1.959964),exp(cad$AME+cad$SE*1.959964),exp(bm$AME),exp(bm$AME-bm$SE*1.959964),exp(bm$AME+bm$SE*1.959964)))

g=ggplot(cad,aes(col=rs1967309_fac,x=endpoint,y=exp(AME)))+geom_point(position = position_dodge(0.7), size=2)+geom_pointrange(size=1.3,aes(x=endpoint,ymin=exp(AME-SE*1.959964),ymax=exp(AME+SE*1.959964)),position = position_dodge(0.7))+geom_hline(yintercept=1,linetype='dashed',color='grey30',size=2)+ coord_flip()+facet_wrap(~male)+colScale+ylim(c(0.75,1.25))


g=g+ theme_bw()+ theme(strip.text.x=element_blank(),panel.spacing=unit(2,'lines'),panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),  panel.border=element_rect(colour='black', linetype='solid'), axis.text.x=element_text(color='black'),axis.text.y=element_text(color='black'),legend.position="none", axis.text=element_text(color='black'))


print(g)

png('UK.BM.png', height=380,width=1200)
xlim=range(c(exp(cad$AME),exp(cad$AME-cad$SE*1.959964),exp(cad$AME+cad$SE*1.959964),exp(bm$AME),exp(bm$AME-bm$SE*1.959964),exp(bm$AME+bm$SE*1.959964)))

g=ggplot(bm,aes(col=rs1967309_fac,x=endpoint,y=exp(AME)))+geom_point(position = position_dodge(0.7), size=2)+geom_pointrange(size=1.3,aes(x=endpoint,ymin=exp(AME-SE*1.959964),ymax=exp(AME+SE*1.959964)),position = position_dodge(0.7))+geom_hline(yintercept=1,linetype='dashed',color='grey30',size=2)+ coord_flip()+facet_wrap(~male)+colScale+ylim(c(0.9,1.1))


g=g+ theme_bw()+ theme(strip.text.x=element_blank(),panel.spacing=unit(2,'lines'),panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),  panel.border=element_rect(colour='black', linetype='solid'), axis.text.x=element_text(color='black'),axis.text.y=element_text(color='black'),legend.position="none", axis.text=element_text(color='black'))


print(g)
