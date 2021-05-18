#!/usr/bin/env Rscript
library("data.table")

print('Figure 2d')
freadf <- function(...) return(as.data.frame(fread(...)))

bysexe=FALSE

# Ignore if bysexe == FALSE
sexe='Female'


 
if(bysexe){
	png(paste(c("NullDistribution.",sexe,".png"),collapse=""),width=500, height=500)
	if(sexe=='Male'){
		# Value extracted from the r2 (you can get it in the file of PEL)
		value=0.348273
		f='SamePairsLimaa.MALE.txt'
	}else{
		# Value extracted from the r2 (you can get it in the file of PEL)
		value=0.00155372
		f='SamePairsLimaa.FEMALE.txt'
	}
	max=0.35

}else{
	png("NullDistribution.png",width=500, height=500)
	# Value extracted from the r2 (you can get it in the file of PEL)
	value=0.0795666
	max=0.17
	f='SamePairsLimaa.All.txt'
}


data_2=read.table(f,header=TRUE)

print('Number of comparison')
print(nrow(data_2))
lower=data_2[,'R.2']>=value

# Calcul the p-value
print('P-value PEL')
print(sum(lower)/length(lower))
print(sum(lower))
print(length(lower))



Saut<-0.005
Fin<-max(range(data_2[,'R.2'],max))


test=as.data.frame(data_2[,'R.2'])
colnames(test)='R2'
test$Pop='PEL'

hist(data_2[,5], breaks=100,main='',cex=2, ylab='',xlab='', col='grey',xlim=c(0,max))

abline(v=quantile(data_2[,5],0.95), lwd=2,col='grey40')
abline(v=value, col='black',lty=2,lwd=2)



