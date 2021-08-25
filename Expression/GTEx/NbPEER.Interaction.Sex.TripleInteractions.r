#!/usr/bin/env Rscript

print('Figure 5 - Figure Supplement 1')


tissu=c('Artery-Tibial','Skin-SunExposed_Lowerleg')
gene='ENSG00000087237.10'

pdf('Sex.3Way.Interaction.pdf',height=7,width=7)

for (organ in tissu) {
	# Plot p-value from the interaction of rs1967309*rs158477 by sex for each PEER factors added
	data=readRDS(paste(c(organ,gene,'rds'),collapse='.'))
	title=c()

	for (i in unique(data$SexColor)) {
		
		for (y in unique(data$BetaSignPch)) {
			temp=data[data$SexColor==i & data$BetaSignPch==y,]

			plot(x=temp$NbPEER,y=-log10(temp$Pvalue), type='l',ylim=(range(c(0,3,-log10(data$Pvalue)))),xlim=(range(c(0,data$NbPEER))), col=i, yaxt='n', xaxt='n',lwd=4,ylab='',xlab='', bty="n")
			par(new=TRUE)


		}
		if(i=='#6e4bf8'){
			title=paste(c('Male : ',unique(temp$Number)),collapse='')
		}else{
			title=paste(c(title,'Female :',unique(temp$Number)),collapse=' ')
		}

	}
	
	plot(x=data$NbPEER,y=-log10(data$Pvalue),main=paste(c(organ,title),collapse='\n'), pch=data$BetaSignPch,cex=2.5,lwd=2,ylab='',xlab='',col=data$SexColor, bg=data$Codingrs1967309Color, ylim=(range(c(0,3,-log10(data$Pvalue)))), yaxt='n', xaxt='n')
	axis(2,las=2, cex.axis=1.5)
	axis(1,xlim=range(data$NbPEER), cex.axis=1.5)
	abline(h=c(-log10(c(0.1,0.05,0.01))),col=c('orange','red','orchid1'), lwd=2, lty=c(2,1,4))


	# Plot p-value from the interaction of rs1967309*rs158477*sex for each PEER factors added
	data2=readRDS(paste(c(organ,gene,'3Way.rds'),collapse='.'))

	print(head(data2))
	plot(x=data2$NbPEER,y=-log10(data2$Pvalue),main=paste(c(organ,sum(unique(data$Number))),collapse='\n'), cex=2, ylab='',xlab='',  ylim=(range(c(0,3,-log10(data2$Pvalue)))), yaxt='n', xaxt='n', pch=18)
	axis(2,las=2, cex.axis=1.5)
	axis(1,xlim=range(data2$NbPEER), cex.axis=1.5)
	abline(h=c(-log10(c(0.1,0.05,0.01))),col=c('orange','red','orchid1'), lwd=2, lty=c(2,1,4))

}