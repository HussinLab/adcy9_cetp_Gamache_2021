#!/usr/bin/env Rscript
library(plyr)
library(reshape2)
library("gplots")
library("ggplot2")
library("RColorBrewer")
library("gridGraphics")
library("grid")
library("ggpubr")
library('MASS')
library('margins')
#Faire scripts dans ../../../Expression/Expression_GTEx/Data/vcf/

pdf('EstimatePheno.pdf', height=2.5,width=6)
table=readRDS('Table.GTEx.rds')

col_7309=c('#56bb5d','#317d36','#0b3c0e')
names(col_7309)=c('GG rs1967309','AG rs1967309','AA rs1967309')
colScale=scale_colour_manual(name='Geno', values=col_7309)

table$Geno=factor(table$Geno,c('GG rs1967309','AG rs1967309','AA rs1967309'))

print(table)
table=table[1:3,]

g=ggplot(table,aes(col=Geno,x=Pheno,y=exp(Estimate)))+geom_point(position = position_dodge(0.7), size=2)+geom_pointrange(size=1.3,aes(x=Pheno,ymin=exp(ymin),ymax=exp(ymax)),position = position_dodge(0.7))+geom_hline(yintercept=1,linetype='dashed',color='grey30',size=2)+ coord_flip()+colScale


g=g+ theme_bw()+ theme(strip.text.x=element_blank(),panel.spacing=unit(2,'lines'),panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),  panel.border=element_rect(colour='black', linetype='solid'), axis.text.x=element_text(color='black'),axis.text.y=element_text(color='black'), legend.position='none', axis.text=element_text(color='black'))
print(g)