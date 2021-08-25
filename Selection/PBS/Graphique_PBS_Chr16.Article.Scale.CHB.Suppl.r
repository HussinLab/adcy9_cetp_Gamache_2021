#!/usr/bin/env Rscript
library("data.table")
library("scales")

print('Figure child 4-1')


freadf <- function(...) return(as.data.frame(fread(...)))

##########################################
#########						##########
#########		Data 			##########
#########						##########
##########################################


diviseur=1000000

ylim=c(-0.02,0.17)
xlim=c(4012650,4166186)

# Conversion Bp to Mb
xlim=xlim/diviseur

# Color used
cPEL="#ff8c65" # #CC79A7
cMXL="#a6d624" # #E69F00
cCHB="#3576ff" # #0072B2


nameY='Branch length of subpopulations'
nameX='Positions on the chromosome 16 (bp)'

main='ADCY9'

exonADCY9=read.table('../iHS/Exon.ADCY9.txt',header=TRUE)
exonADCY9=exonADCY9/diviseur
exonADCY9$ydown=-0.02
exonADCY9$yup=-0.01

# for(s in c('All','MALE','FEMALE')){
for(s in c('MALE','FEMALE')){
	print('Begin read')
	# Coordinate kept from all chromosomes in readed files : chr16; 4002650-4176186
	if(s=='All'){
		dataNormZ=freadf('PBS.ADCY9.txt')
		pdf('PBS.ADCY9.Article.Scale.CHB.pdf',width=11, height=6)

		# Quantile were calculated from all chromosomes,
		# qPEL=quantile(data$PEL,0.95)
		# qMXL=quantile(data$MXL,0.95)
		# qCHB=quantile(data$CHB,0.95)
		qPEL=0.03426227
		qMXL=0.02952179
		qCHB=0.07999562

	}else{

		dataNormZ=freadf(paste(c('PBS.ADCY9.',s,'.txt'),collapse=''))
		# Quantile were calculated from chromosome 16
		# print(quantile(dataNormZ$PEL,0.95))
		# print(quantile(dataNormZ$MXL,0.95))
		# print(quantile(dataNormZ$CHB,0.95))

		if(s=='MALE'){
			pdf('PBS.ADCY9.Article.Scale.CHB.Sex.pdf',width=7.2, height=7.5)
			qPEL=0.04887756
			qMXL=0.04174335 
			qCHB=0.09920801 
		}else{ #Female
			qPEL=0.04543333 
			qMXL=0.03985739  
			qCHB=0.1021854 
		}
	}
	

	colnames(dataNormZ)=c('Position','PEL','MXL','CHB','Chr')


	# To keep the position of the SNP without conversion Bp to Mb
	rownames(dataNormZ)=dataNormZ$Position

	# Conversion Bp to Mb
	dataNormZ$Position=dataNormZ$Position/diviseur


	# Plot all values in ADCY9 gene
	plot(x=dataNormZ$Position,y=dataNormZ$CHB, xlim=xlim, col=cCHB,ylim=ylim, xlab='',ylab='',yaxt='n', lwd=1, bty='n', cex.axis=2, cex.lab=2)
	axis(2,ylim=ylim,las=2, cex.axis=2)
	par(new=TRUE)
	plot(x=dataNormZ$Position,y=dataNormZ$PEL, xlim=xlim, pch=16, col=cPEL,ylim=ylim, xlab='',ylab='', main='',xaxt="n",yaxt="n", lwd=2, bty='n',cex.lab=2, cex=2)
	par(new=TRUE)
	plot(x=dataNormZ$Position,y=dataNormZ$MXL, xlim=xlim, col=cMXL,ylim=ylim, xlab='',ylab='',xaxt="n",yaxt="n", lwd=1, bty='n')

	# Plot value for rs1967309
	# par(new=TRUE)
	# plot(x=4065583/diviseur,y=dataNormZ['4065583','MXL'], xlim=xlim,pch=22, col='black',ylim=ylim, xlab='',ylab='',cex=2,xaxt="n",yaxt="n",lwd=2, bty='n')
	# par(new=TRUE)
	# plot(x=4065583/diviseur,y=dataNormZ['4065583','CHB'], xlim=xlim,pch=23, col='black',ylim=ylim, xlab='',ylab='',cex=2,xaxt="n",yaxt="n",lwd=2, bty='n')
	par(new=TRUE)
	plot(x=4065583/diviseur,y=dataNormZ['4065583','PEL'], xlim=xlim, col='black',ylim=ylim, xlab='',ylab='',cex=2.4,xaxt="n",yaxt="n",lwd=4, bg=cPEL)

	# Added the gene under the graph
	par(new=TRUE)
	segments(4012650/diviseur,-0.015,4166186/diviseur,col='blue') #All genes
	par(new=TRUE)
	rect(exonADCY9$Begin,exonADCY9$ydown,exonADCY9$End,exonADCY9$yup, col='blue',border = 'blue') #Box for exon
	par(new=TRUE)
	rect(4012650/diviseur,-0.013,min(exonADCY9$Begin),-0.017, col='blue',border = 'blue') #UTR
	par(new=TRUE)
	rect(4166186/diviseur,-0.013,max(exonADCY9$End),-0.017, col='blue',border = 'blue') #UTR

	# Add percentile for each population
	abline(h=qPEL,col=cPEL, lwd=2)
	abline(h=qMXL,col=cMXL, lwd=2)
	abline(h=qCHB,col=cCHB, lwd=2)

	# Uncomment if you want the legend
	# legend('topright', c("Peru (PEL)","Mexico (MXL)","East Asia (CHB)    "), col=c(cPEL,cMXL,cCHB), bg="white", pch=c(24,22,23), cex=2, pt.cex=2,pt.lwd=4,bty = "n")

	# Highest summit of the Hotspot around rs1967309
	abline(v=4045116/diviseur, lwd=2)
	abline(v=4077442/diviseur, lwd=2)


	# print()

}