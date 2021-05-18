#!/usr/bin/env Rscript
library("data.table")
library("plyr")
library("scales")
print('Figure 1c')


freadf <- function(...) return(as.data.frame(fread(...)))

##########################################
#########						##########
#########		Data 			##########
#########						##########
##########################################


diviseur=1000000

ylim=c(0,0.17)
	
sexe=FALSE
g2=''
		

chr16=TRUE


print('Begin read')
dataNormZ=freadf('PBS.ADCY9.txt')

colnames(dataNormZ)=c('Position','PEL','MXL','CHB','Chr')

png('PBS.ADCY9.Article.Scale.CHB.png',width=1000, height=500)


# dataNormZ=dataNormZ[dataNormZ$Chr==16 & dataNormZ$Position>=(4012650-10000) & dataNormZ$Position<=(4166186+10000),]
# write.table(dataNormZ,'PBS.ADCY9.txt',col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')

xlim=c(4012650,4166186)

nameY='Branch length of subpopulations'
nameX='Positions on the chromosome 16 (bp)'

main='ADCY9'

par(mar=c(5, 7, 5, 5))

# Bp to Mb
xlim=xlim/diviseur

rownames(dataNormZ)=dataNormZ$Position
dataNormZ$Position=dataNormZ$Position/diviseur

col=c("#56B4E9","#CC79A7","#D55E00","#009E73")

# Plot line by population
plot(x=dataNormZ$Position,y=dataNormZ$CHB, xlim=xlim,type='l', col='#0072B2',ylim=ylim, xlab='',ylab='',yaxt='n', lwd=2, bty='n', cex.axis=2, cex.lab=2)
axis(2,ylim=ylim,las=2, cex.axis=2)
par(new=TRUE)
plot(x=dataNormZ$Position,y=dataNormZ$PEL, xlim=xlim,type='l', col='#CC79A7',ylim=ylim, xlab='',ylab='', main='',xaxt="n",yaxt="n", lwd=2, bty='n',cex.lab=2)
par(new=TRUE)
plot(x=dataNormZ$Position,y=dataNormZ$MXL, xlim=xlim,type='l', col='#E69F00',ylim=ylim, xlab='',ylab='',xaxt="n",yaxt="n", lwd=2, bty='n')

# Plot value for rs1967309
par(new=TRUE)
plot(x=4065583/diviseur,y=dataNormZ['4065583','MXL'], xlim=xlim,pch=22, col='#E69F00',ylim=ylim, xlab='',ylab='',bg='#F0E442',cex=3,xaxt="n",yaxt="n",lwd=3, bty='n')
par(new=TRUE)
plot(x=4065583/diviseur,y=dataNormZ['4065583','CHB'], xlim=xlim,pch=23, col='#0072B2',ylim=ylim, xlab='',ylab='',bg='#F0E442',cex=3,xaxt="n",yaxt="n",lwd=3, bty='n')
par(new=TRUE)
plot(x=4065583/diviseur,y=dataNormZ['4065583','PEL'], xlim=xlim,pch=24, col='#CC79A7',ylim=ylim, xlab='',ylab='',bg='#F0E442',cex=3,xaxt="n",yaxt="n",lwd=3)


# Quantile were calculated from all chromosomes 

# abline(h=quantile(data$PEL,0.95),col='#CC79A7', lwd=2)
# abline(h=quantile(data$MXL,0.95),col='#E69F00', lwd=2)
# abline(h=quantile(data$CHB,0.95),col='#0072B2', lwd=2)

abline(h=0.03426227,col='#CC79A7', lwd=2)
abline(h=0.02952179,col='#E69F00', lwd=2)
abline(h=0.07999562,col='#0072B2', lwd=2)


legend('topright', c("Peru (PEL)","Mexico (MXL)","East Asia (CHB)    "), col=c("#CC79A7","#E69F00","#0072B2"), bg="white", pch=c(24,22,23), cex=2, pt.cex=2,pt.lwd=4,bty = "n")

# Highest summit of the Hotspot around rs1967309
abline(v=4045116/diviseur, lwd=2)
abline(v=4077442/diviseur, lwd=2)




