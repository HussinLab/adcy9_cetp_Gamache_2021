#!/usr/bin/env Rscript
library('scales')
print('Figure 3a,b and 4')


# 1 = rs1967309 in CETP - ALL
# 2 = Top 4 SNPs in ADCY9 - ALL
# 3 = rs1967309 in CETP - Male
# 4 = Top 4 SNPs in ADCY9 - Male
# 5 = rs1967309 in CETP - Female
# 6 = Top 4 SNPs in ADCY9 - Female

plot=4

if(plot==2){
	pdf(paste(c('r2.',plot,'.pdf'),collapse=''),width=15, height=7.5)

}else if(plot%in%c(1,3,5)){
	pdf(paste(c('r2.',plot,'.pdf'),collapse=''),width=7.5, height=7.5)

}else{
	pdf(paste(c('r2.',plot,'.pdf'),collapse=''),width=7.5, height=7)
}

rs1967309=4065583
rs158477=57007610

# From Bp to Mb
diviseur=1000000

if(plot>=3){ 
	par(mar=c(3, 6, 3, 2))

}
xlim=c(4012650,4166186)/diviseur

# All Peruvians
data_all<-read.table("ADCY9_CETP.PEL.005.geno.ld", header=TRUE, row.names=NULL,  sep="\t")
Couleur_6<-c("#56B4E9","#CC79A7","#D55E00","#009E73")

# Get the 99th quantile of all comparisons
percentile_99_All<-quantile(data_all[,5],0.99)
# Extract top SNPs
data_all_CETP=data_all[data_all$POS2 %in% c(57007610,57008227,57008287,57014319),]
# Extract rs1967309
data_all_7309=data_all[data_all$POS1==rs1967309,]

# Assign color
data_all_CETP$Color='white'
data_all_CETP[data_all_CETP$POS2==57007610,'Color']=Couleur_6[1]
data_all_CETP[data_all_CETP$POS2 %in% c(57008227,57008287),'Color']=Couleur_6[2]
data_all_CETP[data_all_CETP$POS2 == 57014319,'Color']=Couleur_6[4]

data_all_7309$Color='black'
data_all_7309[data_all_7309$POS2==57007610,'Color']=Couleur_6[1]
data_all_7309[data_all_7309$POS2 %in% c(57008227,57008287),'Color']=Couleur_6[2]
data_all_7309[data_all_7309$POS2 == 57014319,'Color']=Couleur_6[4]

# Get coordinate from both genes
exonADCY9=read.table('../../iHS/Exon.ADCY9.txt',header=TRUE)

exonADCY9=exonADCY9/diviseur
print(head(exonADCY9))
exonADCY9$ydown=-0.022
exonADCY9$yup=-0.007

exonCETP=read.table('../../Exon.CETP.txt',header=FALSE)
exonCETP=exonCETP[,c(2,3)]
colnames(exonCETP)=c('Begin','End')
exonCETP=exonCETP/diviseur
exonCETP$ydown=-0.022
exonCETP$yup=-0.007


if(plot==1){ #In CETP
	# Bp to Mb
	xlim=c(56995835,57100756)/diviseur
	data_all_7309[,3]=data_all_7309[,3]/diviseur
	plot(x=data_all_7309[,3], y=data_all_7309[,5], pch=16, xlab="CETP", ylab="", main="", col=data_all_7309$Color,ylim=c(-0.02,0.17),yaxt='n',cex=2, cex.axis=2)

	axis(2,ylim=range(data_all_7309[,5]),las=2, cex.axis=2)
	abline(h=percentile_99_All)
	# legend('topright', c('rs158477','rs158480/rs158617','rs12447620'),col=Couleur_6[c(1,2,4)],pch=16,bg="white",bty = "n",cex=1.2)

	par(new=TRUE)
	segments(56995834/diviseur,-0.015,57017756/diviseur,col='blue')
	par(new=TRUE)
	rect(exonCETP$Begin,exonCETP$ydown,exonCETP$End,exonCETP$yup, col='blue',border = 'blue')
	par(new=TRUE)
	rect(56995834/diviseur,-0.013,min(exonCETP$Begin),-0.017, col='blue',border = 'blue')
	par(new=TRUE)
	rect(57017756/diviseur,-0.013,max(exonCETP$End),-0.017, col='blue',border = 'blue')

	print(max(exonCETP$End))



}else if(plot==2){ #In ADCY
	data_all_CETP[,2]=data_all_CETP[,2]/diviseur
	plot(x=data_all_CETP[,2], y=data_all_CETP[,5], pch=16, xlab="ADCY9", ylab="", main="", col=data_all_CETP$Color,ylim=c(-0.02,0.17),yaxt='n',cex=2, cex.axis=2, xlim=xlim)
	axis(2,ylim=range(data_all_CETP[,5]),las=2,cex.axis=2)
	abline(h=percentile_99_All, col='black')

	abline(v=rs1967309/diviseur,lty=2,lwd=2,col='black')
	# legend('topright', c('rs158477','rs158480/rs158617','rs12447620'),col=Couleur_6[c(1,2,4)],pch=16,bg="white",bty = "n",cex=1.2)
	legend('topleft', c(''),lty=2,lwd=2,col='black',cex=1.2,bg="white",bty = "n")

	par(new=TRUE)
	segments(4012650/diviseur,-0.015,4166186/diviseur,col='blue')
	par(new=TRUE)
	rect(exonADCY9$Begin,exonADCY9$ydown,exonADCY9$End,exonADCY9$yup, col='blue',border = 'blue')
	par(new=TRUE)
	rect(4012650/diviseur,-0.013,min(exonADCY9$Begin),-0.017, col='blue',border = 'blue')
	par(new=TRUE)
	rect(4166186/diviseur,-0.013,max(exonADCY9$End),-0.017, col='blue',border = 'blue')

}








#Male Peruvians
data_men<-read.table("ADCY9_CETP.PEL.005.MALE.geno.ld", header=TRUE, row.names=NULL,  sep="\t")
data_men_CETP=data_men[data_men$POS2 %in% c(57007610,57008227,57008287,57014319),]
data_men_7309=data_men[data_men$POS1==rs1967309,]
percentile_99_men<-quantile(data_men[,5],0.99)



data_men_CETP$Color=NA
data_men_CETP[data_men_CETP$POS2==57007610,'Color']=Couleur_6[1]
data_men_CETP[data_men_CETP$POS2 %in% c(57008227),'Color']=Couleur_6[2]
data_men_CETP[data_men_CETP$POS2 %in% c(57008287),'Color']=Couleur_6[2]
data_men_CETP[data_men_CETP$POS2 == 57014319,'Color']=Couleur_6[2]


data_men_7309$Color='black'
data_men_7309[data_men_7309$POS2==57007610,'Color']=Couleur_6[1]
data_men_7309[data_men_7309$POS2 %in% c(57008227,57008287),'Color']=Couleur_6[2]
data_men_7309[data_men_7309$POS2 == 57014319,'Color']=Couleur_6[4]

pos=sort(data_men_7309$POS2)
pos=c(pos[!(pos%in%c(57007610,57008227,57008287,57014319))],57007610,57008227,57008287,57014319)

data_men_7309=data_men_7309[match(pos,data_men_7309$POS2),]


if(plot==3){
	data_men_7309[,3]=data_men_7309[,3]/diviseur
	plot(x=data_men_7309[,3], y=data_men_7309[,5], pch=16, xlab="CETP", ylab="r2 value", main="", col=data_men_7309$Color,ylim=c(-0.02,0.37), yaxt='n',cex=1.2)
	axis(2,ylim=range(data_men_7309[,5]),las=1)
	abline(h=percentile_99_men<-quantile(data_men[,5],0.99))
	legend('topright', c('rs158477','rs158480/rs158617','rs12447620'),col=c(Couleur_6[c(1,2,4)]),pch=16,bg="white",bty = "n",cex=1.5)

	par(new=TRUE)
	segments(56995834/diviseur,-0.015,57017756/diviseur,col='blue')
	par(new=TRUE)
	rect(exonCETP$Begin,exonCETP$ydown,exonCETP$End,exonCETP$yup, col='blue',border = 'blue')
	par(new=TRUE)
	rect(56995834/diviseur,-0.013,min(exonCETP$Begin),-0.017, col='blue',border = 'blue')
	par(new=TRUE)
	rect(57017756/diviseur,-0.013,max(exonCETP$End),-0.017, col='blue',border = 'blue')
	

}else if(plot==4){

	data_men_CETP[,2]=data_men_CETP[,2]/diviseur
	plot(x=data_men_CETP[,2], y=data_men_CETP[,5], pch=16, xlab="ADCY9", ylab="", main="", col=data_men_CETP$Color,ylim=c(-0.02,0.37),yaxt='n',cex=2, cex.axis=2, xlim=xlim)
	axis(2,ylim=range(data_men_CETP[,5]),las=2,cex.axis=2)
	abline(h=percentile_99_men)
	# legend('topright', c('rs158477','rs158480/rs158617/rs12447620'),col=c(Couleur_6[c(1,2)]),pch=16,bg="white",bty = "n",cex=1.2)
	legend('topleft', c(''),lty=2,lwd=2,bg="white",bty = "n",cex=1.2)
	abline(v=rs1967309/diviseur,lty=2,lwd=2,col='black')

	par(new=TRUE)
	segments(4012650/diviseur,-0.015,4166186/diviseur,col='blue')
	par(new=TRUE)
	rect(exonADCY9$Begin,exonADCY9$ydown,exonADCY9$End,exonADCY9$yup, col='blue',border = 'blue')
	par(new=TRUE)
	rect(4012650/diviseur,-0.013,min(exonADCY9$Begin),-0.017, col='blue',border = 'blue')
	par(new=TRUE)
	rect(4166186/diviseur,-0.013,max(exonADCY9$End),-0.017, col='blue',border = 'blue')

}





# Female Peruvian
data_men<-read.table("ADCY9_CETP.PEL.005.FEMALE.geno.ld", header=TRUE, row.names=NULL,  sep="\t")
data_men_CETP=data_men[data_men$POS2 %in% c(57007610,57008227,57008287,57014319),]
data_men_7309=data_men[data_men$POS1==rs1967309,]
percentile_99_men<-quantile(data_men[,5],0.99)


data_men_CETP$Color=NA
data_men_CETP[data_men_CETP$POS2==57007610,'Color']=Couleur_6[1]
data_men_CETP[data_men_CETP$POS2 %in% c(57008227,57008287),'Color']=Couleur_6[2]
data_men_CETP[data_men_CETP$POS2 == 57014319,'Color']=Couleur_6[2]


data_men_7309$Color='black'
data_men_7309[data_men_7309$POS2==57007610,'Color']=Couleur_6[1]
data_men_7309[data_men_7309$POS2 %in% c(57008227,57008287),'Color']=Couleur_6[2]
data_men_7309[data_men_7309$POS2 == 57014319,'Color']=Couleur_6[4]

pos=sort(data_men_7309$POS2)
pos=c(pos[!(pos%in%c(57007610,57008227,57008287,57014319))],57007610,57008227,57008287,57014319)

data_men_7309=data_men_7309[match(pos,data_men_7309$POS2),]

if(plot==5){
	data_men_7309[,3]=data_men_7309[,3]/diviseur
	plot(x=data_men_7309[,3], y=data_men_7309[,5], pch=16, xlab="Position in CETP (bp)", ylab="r2 value", main="Peru (PEL)\nWomen", col=data_men_7309$Color,ylim=c(-0.02,0.37), yaxt='n',cex=1.2)
	#legend('topright', 'E', box.col=NA, cex=1.2,bg="white")
	axis(2,ylim=range(data_men_7309[,5]),las=1)
	abline(h=percentile_99_men<-quantile(data_men[,5],0.99))

	legend('topright', c('rs158477','rs158480/rs158617','rs12447620'),col=Couleur_6[c(1,2,4)],pch=16,bg="white",bty = "n",cex=1.5)

	par(new=TRUE)
	segments(56995834/diviseur,-0.015,57017756/diviseur,col='blue')
	par(new=TRUE)
	rect(exonCETP$Begin,exonCETP$ydown,exonCETP$End,exonCETP$yup, col='blue',border = 'blue')
	par(new=TRUE)
	rect(56995834/diviseur,-0.013,min(exonCETP$Begin),-0.017, col='blue',border = 'blue')
	par(new=TRUE)
	rect(57017756/diviseur,-0.013,max(exonCETP$End),-0.017, col='blue',border = 'blue')


}else if(plot==6){

	data_men_CETP[,2]=data_men_CETP[,2]/diviseur
	plot(x=data_men_CETP[,2], y=data_men_CETP[,5], pch=16, xlab="ADCY9", ylab="", main="", col=data_men_CETP$Color,ylim=c(-0.02,0.37),yaxt='n',cex=2,cex.axis=2, xlim=xlim)
	axis(2,ylim=range(data_men_CETP[,5]),las=2,cex.axis=2)
	abline(h=percentile_99_men)
	abline(v=rs1967309/diviseur,lty=2,lwd=2,col='black')

	# legend('topright', c('rs158477','rs158480/rs158617/rs12447620'),col=c(Couleur_6[c(1,2)]),pch=16,bg="white",bty = "n",cex=1.2)
	legend('topleft', c(''),lty=2,lwd=2,bg="white",bty = "n",cex=1.2)

	par(new=TRUE)
	segments(4012650/diviseur,-0.015,4166186/diviseur,col='blue')
	par(new=TRUE)
	rect(exonADCY9$Begin,exonADCY9$ydown,exonADCY9$End,exonADCY9$yup, col='blue',border = 'blue')
	par(new=TRUE)
	rect(4012650/diviseur,-0.013,min(exonADCY9$Begin),-0.017, col='blue',border = 'blue')
	par(new=TRUE)
	rect(4166186/diviseur,-0.013,max(exonADCY9$End),-0.017, col='blue',border = 'blue')


}

dev.off()