#!/usr/bin/env Rscript
library('scales')
print('Figure 4 (and supplement figures)')

# Position of my main SNP
rs1967309=4065583

# Conversion Bp to Mb 
diviseur=1000000

# each sex
subpop=c('all','Men','Women')

# ADCY9 and CETP genes
exonADCY9=read.table('../../iHS/Exon.ADCY9.txt',header=TRUE) #Exon coordinate from UCSC
exonADCY9=exonADCY9/diviseur

exonCETP=read.table('../../Exon.CETP.txt',header=FALSE) #Raw file get from UCSC
exonCETP=exonCETP[,c(2,3)]
colnames(exonCETP)=c('Begin','End') 
exonCETP=exonCETP/diviseur


# To draw the gene under
exonADCY9$ydown=-0.0004
exonADCY9$yup=-0.001

exonCETP$ydown=-0.0004
exonCETP$yup=-0.001

# Color for SNPs higher than the 99th percentile for LRLD in PEL
Couleur_6<-c("#56B4E9","#CC79A7","#D55E00","#009E73")

for (plot in seq(1, 6, 1)) {

	pdf(paste(c("r2.LIMAA.",plot,".pdf"),collapse=''),width=7.5, height=7)
	par(mar=c(3, 6, 3, 2))


	ylim=c(-0.001,0.01)
	pop=subpop[1]



	if(plot==3 || plot==4){
		pop=subpop[2]
		ylim=c(-0.001,0.018)
	}

	if(plot==5 || plot==6){
		pop=subpop[3]
		ylim=c(-0.001,0.018)
	}

	

	if(pop=='all'){
		data_all<-read.table(paste(c("All.005.BigBlob.geno.ld"),collapse=''), header=TRUE, row.names=NULL,  sep="\t")
	}else{
		data_all<-read.table(paste(c("All.",pop,".005.BigBlob.geno.ld"),collapse=''), header=TRUE, row.names=NULL,  sep="\t")
	}
	
	

	# Get the 99th percentile of all comparisons
	percentile_99_All<-quantile(data_all[,5],0.99)
	
	data_all_7309=data_all[data_all$POS1==rs1967309,]
	data_CETP=data_all[data_all$POS2 %in% c(57007610,57008227,57008287,57014319),]


	# extract data for rs1967309 and color for each SNPs above 99th percentile of LRLD
	data_all_7309$Color='black' #None should be black at the end
	data_all_7309[data_all_7309$POS2==57007610,'Color']=Couleur_6[1]
	data_all_7309[data_all_7309$POS2 %in% c(57008227),'Color']=Couleur_6[2]
	data_all_7309[data_all_7309$POS2 %in% c(57008287),'Color']=Couleur_6[3]
	data_all_7309[data_all_7309$POS2 == 57014319,'Color']=Couleur_6[4]

	
	# extract data for CETP and color for each SNPs above 99th percentile of LRLD
	data_CETP$Color='black'
	data_CETP[data_CETP$POS2==57007610,'Color']=Couleur_6[1]
	data_CETP[data_CETP$POS2 %in% c(57008227),'Color']=Couleur_6[2]
	data_CETP[data_CETP$POS2 %in% c(57008287),'Color']=Couleur_6[2] #Same color that previous because perfect LD
	data_CETP[data_CETP$POS2 == 57014319,'Color']=Couleur_6[2]

	
	if(plot==1 || plot==3 || plot==5){ # CETP gene

		xlim=c(56995835,57017756)/diviseur

		# So that rs158477 can be in the front
		data_all_477=data_all_7309[data_all_7309$POS2==57007610,]
		data_all_7309=data_all_7309[data_all_7309$POS2!=57007610,]
		data_all_7309=rbind(data_all_7309,data_all_477)

		# Convert Bp to Mb
		data_all_7309[,3]=data_all_7309[,3]/diviseur

		plot(x=data_all_7309[,3], y=data_all_7309[,5], pch=16, xlab="", ylab="", main='', col=data_all_7309$Color,ylim=ylim,yaxt='n',cex=1.2, cex.axis=2,xlim=xlim)
		axis(2,ylim=range(data_all_7309[,5]),las=2, cex.axis=2)
		abline(h=percentile_99_All)

		# Legend for each SNPs in CETP
		# legend('topright', c('rs158477','rs158480','rs158617','rs12447620'),col=Couleur_6,pch=16,bg="white",bty = "n",cex=1.4)

		# Add CETP gene under
		par(new=TRUE)
		segments(56995834/diviseur,-0.0007,57017756/diviseur,col='blue')
		par(new=TRUE)
		rect(exonCETP$Begin,exonCETP$ydown,exonCETP$End,exonCETP$yup, col='blue',border = 'blue')
		par(new=TRUE)
		rect(56995834/diviseur,-0.0008,min(exonCETP$Begin),-0.0006, col='blue',border = 'blue')
		par(new=TRUE)
		rect(57017756/diviseur,-0.0008,max(exonCETP$End),-0.0006, col='blue',border = 'blue')

	}else{ # ADCY9

		# Convert Bp to Mb
		data_CETP[,2]=data_CETP[,2]/diviseur


		plot(x=data_CETP[,2], y=data_CETP[,5], pch=16, xlab="", ylab="", main='', col=data_CETP$Color,ylim=ylim,yaxt='n',cex=2, cex.axis=2)
		axis(2,ylim=range(data_CETP[,5]),las=2, cex.axis=2)
		abline(h=percentile_99_All)
		abline(v=rs1967309/diviseur,lty=2,lwd=2,col='black')

		# Legend for each SNPs in CETP
		# legend('topright', c('rs158477','rs158480/rs158617/rs12447620'),col=c(Couleur_6[c(1,2)]),pch=16,bg="white",bty = "n",cex=1.2)
		legend('topleft', c(''),lty=2,lwd=2,bg="white",bty = "n",cex=1.2) #rs1967309

		# Add ADCY9 gene under
		par(new=TRUE)
		segments(4012650/diviseur,-0.0007,4166186/diviseur,col='blue')
		par(new=TRUE)
		rect(exonADCY9$Begin,exonADCY9$ydown,exonADCY9$End,exonADCY9$yup, col='blue',border = 'blue')
		par(new=TRUE)
		rect(4012650/diviseur,-0.0008,min(exonADCY9$Begin),-0.0006, col='blue',border = 'blue')
		par(new=TRUE)
		rect(4166186/diviseur,-0.0008,max(exonADCY9$End),-0.0006, col='blue',border = 'blue')

		
	}

	dev.off()
}