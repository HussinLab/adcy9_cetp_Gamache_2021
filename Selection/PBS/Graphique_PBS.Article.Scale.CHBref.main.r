#!/usr/bin/env Rscript
library("data.table")
library("scales")

print('Figure 2c & 4a,b')


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


print('Begin read')
# Coordinate kept from all chromosomes in readed files : chr16; 4002650-4176186

dataNormZ=freadf('PBS.ADCY9.txt')
pdf('PBS.ADCY9.Article.Scale.CHBref.SplitPop.pdf',width=11, height=5.5)

par(mfrow=c(1,3))

# Quantile were calculated from all chromosomes,
# qPEL=quantile(data$PEL,0.95)
# qMXL=quantile(data$MXL,0.95)
# qCHB=quantile(data$CHB,0.95)
qPEL=0.03426227
qMXL=0.02952179
qCHB=0.07999562



par(mar = c(2, 0.3, 2, 0.5))
colnames(dataNormZ)=c('Position','PEL','MXL','CHB','Chr')

# To keep the position of the SNP without conversion Bp to Mb
rownames(dataNormZ)=dataNormZ$Position

# Conversion Bp to Mb
dataNormZ$Position=dataNormZ$Position/diviseur


# Plot all values in ADCY9 gene

# Subgraph for the CHB population

# the 0.1 is the max value that we permit, but CHB can go higher
# If you want all values for a population, ylab will differs
# This subscript will uniformise the size of the ADCY9 gene under the graph
m=-0.1/20
ylim=range(c(m-m/30,0.1))
s=0.1*0.03
exonADCY9$ydown=m-s
exonADCY9$yup=m+s

plot(x=dataNormZ$Position,y=dataNormZ$CHB, xlim=xlim, pch=16,col=cCHB,ylim=ylim, xlab='',ylab='',yaxt='n', lwd=2, bty='n', cex.axis=2, cex.lab=2,cex=2)
axis(2,las=2, cex.axis=2)
par(new=TRUE)
plot(x=4065583/diviseur,y=dataNormZ['4065583','CHB'], xlim=xlim, col='black',ylim=ylim, xlab='',ylab='',cex=2.4,xaxt="n",yaxt="n",lwd=4, bg=cCHB)
abline(h=qCHB,col=cCHB, lwd=2)

par(new=TRUE)
segments(4012650/diviseur,m,4166186/diviseur,col='blue') #All genes
par(new=TRUE)
rect(exonADCY9$Begin,exonADCY9$ydown,exonADCY9$End,exonADCY9$yup, col='blue',border = 'blue') #Box for exon
par(new=TRUE)
rect(4012650/diviseur,m+s/2,min(exonADCY9$Begin),m-s/2, col='blue',border = 'blue') #UTR
par(new=TRUE)
rect(4166186/diviseur,m+s/2,max(exonADCY9$End),m-s/2, col='blue',border = 'blue') #UTR
abline(v=4045116/diviseur, lwd=2, lty=2)
abline(v=4077442/diviseur, lwd=2, lty=2)



# Subgraph for the PEL population
m=-0.1/20
ylim=range(c(m-m/30,dataNormZ$PEL,0.1))
s=0.1*0.03
exonADCY9$ydown=m-s
exonADCY9$yup=m+s
plot(x=dataNormZ$Position,y=dataNormZ$PEL, xlim=xlim, pch=16, col=cPEL,ylim=ylim, xlab='',ylab='', main='', lwd=2, bty='n',yaxt="n",cex.axis=2, cex.lab=2, cex=2)
par(new=TRUE)
plot(x=4065583/diviseur,y=dataNormZ['4065583','PEL'], xlim=xlim, col='black',ylim=ylim, xlab='',ylab='',cex=2.4,lwd=4, xaxt="n",yaxt="n",bg=cPEL)
abline(h=qPEL,col=cPEL, lwd=2)

par(new=TRUE)
segments(4012650/diviseur,m,4166186/diviseur,col='blue') #All genes
par(new=TRUE)
rect(exonADCY9$Begin,exonADCY9$ydown,exonADCY9$End,exonADCY9$yup, col='blue',border = 'blue') #Box for exon
par(new=TRUE)
rect(4012650/diviseur,m+s/2,min(exonADCY9$Begin),m-s/2, col='blue',border = 'blue') #UTR
par(new=TRUE)
rect(4166186/diviseur,m+s/2,max(exonADCY9$End),m-s/2, col='blue',border = 'blue') #UTR

abline(v=4045116/diviseur, lwd=2, lty=2)
abline(v=4077442/diviseur, lwd=2, lty=2)



# Subgraph for the MXL population
m=-0.1/20
ylim=range(c(m-m/30,dataNormZ$MXL,0.1))
s=0.1*0.03
exonADCY9$ydown=m-s
exonADCY9$yup=m+s
plot(x=dataNormZ$Position,y=dataNormZ$MXL, xlim=xlim, pch=16,col=cMXL,ylim=ylim, xlab='',ylab='', lwd=2, yaxt="n",bty='n',cex=2, cex.axis=2, cex.lab=2)
par(new=TRUE)
plot(x=4065583/diviseur,y=dataNormZ['4065583','MXL'], xlim=xlim, col='black',ylim=ylim, xlab='',ylab='',cex=2.4,lwd=4, xaxt="n",yaxt="n",bg=cMXL)
abline(h=qMXL,col=cMXL, lwd=2)
par(new=TRUE)
segments(4012650/diviseur,m,4166186/diviseur,col='blue') #All genes
par(new=TRUE)
rect(exonADCY9$Begin,exonADCY9$ydown,exonADCY9$End,exonADCY9$yup, col='blue',border = 'blue') #Box for exon
par(new=TRUE)
rect(4012650/diviseur,m+s/2,min(exonADCY9$Begin),m-s/2, col='blue',border = 'blue') #UTR
par(new=TRUE)
rect(4166186/diviseur,m+s/2,max(exonADCY9$End),m-s/2, col='blue',border = 'blue') #UTR


# Uncomment if you want the legend
# legend('topright', c("Peru (PEL)","Mexico (MXL)","East Asia (CHB)    "), col=c(cPEL,cMXL,cCHB), bg="white", pch=c(24,22,23), cex=2, pt.cex=2,pt.lwd=4,bty = "n")

# Highest summit of the Hotspot around rs1967309
abline(v=4045116/diviseur, lwd=2, lty=2)
abline(v=4077442/diviseur, lwd=2, lty=2)


	