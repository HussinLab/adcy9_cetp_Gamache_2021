#!/usr/bin/env Rscript
library("data.table")

print('Figure 3d and some never presented by figures (by sex)')
freadf <- function(...) return(as.data.frame(fread(...)))

bysexe=FALSE

# Ignore if bysexe = FALSE
# Can be Male of Female if bysexe=TRUE
sexe='Female'


 
if(bysexe){
	pdf(paste(c("NullDistribution.",sexe,".pdf"),collapse=""),width=7.5, height=7.5)

	if(sexe=='Male'){
		# Value extracted from the r2 (you can get it in the file in ../PEL)
		value=0.348273
		f='SamePairsLimaa.MALE.txt'
	}else{
		# Value extracted from the r2 (you can get it in the file in ../PEL)
		value=0.00155372
		f='SamePairsLimaa.FEMALE.txt'
	}
	max=0.35

}else{
	pdf("NullDistribution.pdf",width=7.5, height=7.5)

	# Value extracted from the r2 (you can get it in the file in ../PEL)
	value=0.0795666
	max=0.17
	f='SamePairsLimaa.All.txt'
}


data=read.table(f,header=TRUE)


print('Number of comparison')
print(nrow(data))


# Calcul and print the p-value
lower=data[,'R.2']>=value
print('P-value PEL')
print(sum(lower)/length(lower))
# print(sum(lower))
# print(length(lower))



hist(data[,5], breaks=100,main='',cex=2, ylab='',xlab='', col='grey',xlim=c(0,max))

abline(v=quantile(data[,5],0.95), lwd=2,col='grey40')
abline(v=value, col='black',lty=2,lwd=2)



