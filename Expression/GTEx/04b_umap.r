#!/usr/bin/env Rscript
library("data.table")
library("plyr")
library('ggplot2')
library('reshape2')
library(scales)
library('umap')

###################################################
####                                          #####
####          Fonctions				          #####
####                                          #####
###################################################

freadf <- function(...) return(as.data.frame(fread(...)))


###################################################
####                                          #####
####          data 					          #####
####                                          #####
###################################################

# pcs obtained from 04a script
data1=freadf('pcs.txt')
pca=data1


###################################################
####                                          #####
####          plots 				          #####
####                                          #####
###################################################



if(!file.exists('GTEx.umap.txt')){
  rownames(data1)=data1[,2]
  data1=data1[,-c(1,2)]
  print('Compute umap - Begin')
  a=umap(data1)$layout
  print('Compute umap - Finish')
  colnames(a)=c('umap1','umap2')
  id=as.data.frame(rownames(a))
  colnames(id)='IID'
  a=cbind(id,a)

  write.table(a,'GTEx.umap.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')

  
}else{
  data1=data1[,-1]
  colnames(data1)[1]='SUBJID'
  a=read.table('GTEx.umap.txt',header=TRUE)

}

age=read.table('../Data/SpecificPhenotype/AGE.txt', header=TRUE)
sex=read.table('../Data/SpecificPhenotype/SEX.txt', header=TRUE)
race=read.table('../Data/SpecificPhenotype/RACE.txt', header=TRUE)
ethncty=read.table('../Data/SpecificPhenotype/ETHNCTY.txt', header=TRUE)


cov=merge(age,sex,by='SUBJID')
cov=merge(cov,race,by='SUBJID')
cov=merge(cov,ethncty,by='SUBJID')

colnames(a)[1]='SUBJID'
a=merge(a,data1,by='SUBJID')


cov$Race='Other/Unknown'
r1=cov$RACE==1
cov[r1,'Race']='Asian'
r1=cov$RACE==2
cov[r1,'Race']='Black/African'
r1=cov$RACE==3
cov[r1,'Race']='White'
r1=cov$ETHNCTY==1
cov[r1,'Race']='Latino'

cov$Race=factor(cov$Race,levels=c('White','Black/African','Asian','Latino','Other/Unknown'))



allData=merge(a,cov,by='SUBJID')


ggplot(allData,aes(x=umap1,y=umap2,color=Race))+geom_point()+scale_color_brewer(palette='Set1')+theme(panel.background=element_rect(fill=NA), panel.border=element_rect(colour='black', fill=NA), axis.text=element_text(size=15, color='black'), legend.position='none')


white=allData$RACE==3
allData=allData[white,]


print('Manual extraction of the white non latino, so value will change here depending on your umap')
white_filter=allData$umap1>-8 & allData$umap2<2.5 
allData=allData[white_filter,]

print(sum(white_filter))

keep=allData$SUBJID
write.table(keep, 'WhiteToKeep.txt',row.names=FALSE,quote=FALSE,col.names=FALSE)

dev.off()
