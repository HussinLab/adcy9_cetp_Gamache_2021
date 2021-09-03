#!/usr/bin/env Rscript
library("data.table")
library('ggplot2')
library('reshape2') 

print('Figure 3-figure supplement 2')

freadf <- function(...) return(as.data.frame(fread(...)))



pdf('Comparaison.pdf',width=10, height=10)


commun=readRDS('Freq.Cohorts.rds')

#Keep position with maf of at least 5% in both cohort
positionToKeep=commun[commun$maf1000G>0.05 & commun$mafLIMAA>0.05,'POS']



r2_1000G=read.table('R2Values.PEL-LIMAA.txt',header=TRUE)

Couleur_6<-c("#56B4E9","#CC79A7","#D55E00","#009E73")


ylim=range(c(r2_1000G$R2PEL,r2_1000G$R2LIMAA))
r2_1000G$R2PEL=as.numeric(r2_1000G$R2PEL)
r2_1000G$R2LIMAA=as.numeric(r2_1000G$R2LIMAA)

r2_1000G=r2_1000G[!is.na(r2_1000G[,1]),]


# Color and put our SNP at the front so that they can be visible
r2_1000G$Color='black'
ourSNPs=(r2_1000G$POS1== 4065583 & r2_1000G$POS2== 57007610) 
r2_1000G[ourSNPs,'Color']=Couleur_6[1]
r2_1000G=rbind(r2_1000G,r2_1000G[ourSNPs,])
r2_1000G=r2_1000G[!ourSNPs,]

ourSNPs=(r2_1000G$POS1== 4065583 & r2_1000G$POS2== 57008227) 
r2_1000G[ourSNPs,'Color']=Couleur_6[2]
r2_1000G=rbind(r2_1000G,r2_1000G[ourSNPs,])
r2_1000G=r2_1000G[!ourSNPs,]

ourSNPs=(r2_1000G$POS1== 4065583 & r2_1000G$POS2== 57008287) 
r2_1000G[ourSNPs,'Color']=Couleur_6[3]
r2_1000G=rbind(r2_1000G,r2_1000G[ourSNPs,])
r2_1000G=r2_1000G[!ourSNPs,]

ourSNPs=(r2_1000G$POS1== 4065583 & r2_1000G$POS2== 57014319) 
r2_1000G[ourSNPs,'Color']=Couleur_6[4]
r2_1000G=rbind(r2_1000G,r2_1000G[ourSNPs,])
r2_1000G=r2_1000G[!ourSNPs,]

# Filter frquence
r2_1000G=r2_1000G[(r2_1000G$POS1 %in% positionToKeep) & (r2_1000G$POS2 %in% positionToKeep), ]

plot(x=r2_1000G$R2PEL,y=r2_1000G$R2LIMAA,pch=16, col=r2_1000G$Color, cex=1.5)

abline(h=quantile(r2_1000G$R2LIMAA,0.99))
abline(v=quantile(r2_1000G$R2PEL,0.99))
legend('topright', c('rs158477','rs158480','rs158617','rs12447620'),col=Couleur_6,pch=16,bg="white",bty = "n",cex=1.5)
	