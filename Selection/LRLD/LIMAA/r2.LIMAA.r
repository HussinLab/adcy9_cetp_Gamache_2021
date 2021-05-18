#!/usr/bin/env Rscript
library('scales')
print('Figure 3 (and more)')


rs1967309=4065583
rs158477=57007610

subpop=c('all','Men','Women')

for (plot in seq(1, 6, 1)) {

	png(paste(c("r2.LIMAA.",plot,".png"),collapse=''),width=550, height=500)
	par(mar=c(3, 6, 3, 2))


	ylim=c(0,0.01)
	pop=subpop[1]
	xlim=c(4012650,4166023)#4166186)

	if(plot==3 || plot==4){
		pop=subpop[2]
		ylim=c(0,0.015)
	}

	if(plot==5 || plot==6){
		pop=subpop[3]
		ylim=c(0,0.015)
	}

	

	if(pop=='all'){
		data_all<-read.table(paste(c("All.005.geno.ld"),collapse=''), header=TRUE, row.names=NULL,  sep="\t")
	}else{
		data_all<-read.table(paste(c("All.",pop,".005.geno.ld"),collapse=''), header=TRUE, row.names=NULL,  sep="\t")
	}
	
	Couleur_6<-c("#56B4E9","#CC79A7","#D55E00","#009E73")

	# Get the 99th percentile of all comparisons
	percentile_99_All<-quantile(data_all[,5],0.99)
	
	data_all_7309=data_all[data_all$POS1==rs1967309,]
	data_CETP=data_all[data_all$POS2 %in% c(57007610,57008227,57008287,57014319),]


	# extract data for rs1967309
	data_all_7309$Color='black'
	data_all_7309[data_all_7309$POS2==57007610,'Color']=Couleur_6[1]
	data_all_7309[data_all_7309$POS2 %in% c(57008227),'Color']=Couleur_6[2]
	data_all_7309[data_all_7309$POS2 %in% c(57008287),'Color']=Couleur_6[3]
	data_all_7309[data_all_7309$POS2 == 57014319,'Color']=Couleur_6[4]

	
	# extract data for CETP
	data_CETP$Color='black'
	data_CETP[data_CETP$POS2==57007610,'Color']=Couleur_6[1]
	data_CETP[data_CETP$POS2 %in% c(57008227),'Color']=Couleur_6[2]
	data_CETP[data_CETP$POS2 %in% c(57008287),'Color']=Couleur_6[2]
	data_CETP[data_CETP$POS2 == 57014319,'Color']=Couleur_6[2]

	
	if(plot==1 || plot==3 || plot==5){
		data_all_477=data_all_7309[data_all_7309$POS2==57007610,]
		data_all_7309=data_all_7309[data_all_7309$POS2!=57007610,]
		data_all_7309=rbind(data_all_7309,data_all_477)
		data_all_7309[,3]=data_all_7309[,3]/1000000
		plot(x=data_all_7309[,3], y=data_all_7309[,5], pch=16, xlab="", ylab="", main='', col=data_all_7309$Color,ylim=ylim,yaxt='n',cex=1.2, cex.axis=2)
		axis(2,ylim=range(data_all_7309[,5]),las=2, cex.axis=2)
		abline(h=percentile_99_All)
		legend('topright', c('rs158477','rs158480','rs158617','rs12447620'),col=Couleur_6,pch=16,bg="white",bty = "n",cex=1.4)

	}else{
		data_CETP[,2]=data_CETP[,2]/1000000
		plot(x=data_CETP[,2], y=data_CETP[,5], pch=16, xlab="", ylab="", main='', col=data_CETP$Color,ylim=ylim,yaxt='n',cex=2, cex.axis=2)
		axis(2,ylim=range(data_CETP[,5]),las=2, cex.axis=2)
		abline(h=percentile_99_All)
		abline(v=rs1967309/1000000,lty=2,lwd=2,col='black')
		legend('topright', c('rs158477','rs158480/rs158617/rs12447620'),col=c(Couleur_6[c(1,2)]),pch=16,bg="white",bty = "n",cex=1.2)
		legend('topleft', c('rs1967309'),lty=2,lwd=2,bg="white",bty = "n",cex=1.2)

		
	}

	dev.off()
}