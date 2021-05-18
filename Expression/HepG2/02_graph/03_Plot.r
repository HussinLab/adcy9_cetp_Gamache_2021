#!/usr/bin/env Rscript
library("data.table")
library("plyr")
library('ggpubr')
library(ggplot2)

print('Figure 4a')

##########################################
######							##########
######		Fonction			##########
######							##########
##########################################


format_pval<-function(x){
	return(formatC(format="e",x,digits=2))
}


for_hepG2KO=function(data){
	i=1
	number=c()
	type=c()

	#Uniform name
	while(i<=ncol(data)){
		t=unlist(strsplit(colnames(data)[i],'.',fixed=TRUE))
		number=c(number,t[3])
		l=length(t)
		type=c(type,paste(c(t[l-1],t[l]),collapse="-"))
		i=i+1
	}

	data=t(data)


	a=as.data.frame(type)
	data=cbind(a,data)
	cout=count(data$type)


	colnames(data)=c('type','ADCY9','CETP')

	data$experiment=c(rep('scr',5),rep('si-1039',5))


	#ADCY9
	print(format_pval(t.test(as.numeric(data$ADCY9)~data$experiment,paired=TRUE)$p.value))
	# CETP
	print(format_pval(t.test(as.numeric(data$CETP)~data$experiment,paired=TRUE)$p.value))

	data$experiment=factor(data$experiment,c('scr','si-1039'))
	
	data$Gene=NA
	data$Expression=NA

	data2=data

	data$Gene='ADCY9'
	data$Expression=data$ADCY9

	data2$Gene='CETP'
	data2$Expression=data2$CETP

	data=rbind(data,data2)


	print(ggplot(data,aes(y=Expression,fill=experiment,x=Gene))+geom_boxplot(outlier.shape=NA)+geom_jitter(position=position_dodge(0.75), size=2)+theme(axis.text=element_text(size=20), panel.background=element_rect(fill='white', color='black'), legend.position=c(0.85,0.92),legend.title=element_blank(),legend.text=element_text(size=15))+scale_fill_manual(values=c('#56B4E9','#D55E00')))
	


}

##########################################
######							##########
######		Variable			##########
######							##########
##########################################


png('hepG2KO.no205.CETP.RSEM.Hg38.png',width=500, height=500)
data_all=read.table(paste(c("genesExpressionTable.Hg38.gene6_5.normalizedVoom.ADCY9_CETP.ScrSi.txt"),collapse=""),row.names=1,header=TRUE)
for_hepG2KO(data_all)



	


	





