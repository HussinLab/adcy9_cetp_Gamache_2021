#!/usr/bin/env Rscript
library('reshape2')
library("data.table")
library('scales')

print('Figure 3/4-figure supplement 1/3a,b')
###################################################
####                                          #####
####          Variable				          #####
####                                          #####
###################################################

freadf <- function(...) return(as.data.frame(fread(..., fill=TRUE)))

diviseur=1000000
subpop=c('andean')

for(plot in seq(1, 6, 1)){

	print(plot)
	pdf(paste(c('r2.NAGD.',plot,'.pdf'),collapse=''),width=9, height=7.5)
	par(mar=c(3, 6, 3, 2))


	rs1967309=4065583
	rs158477=57007610


	xlim=c(4012650,4166186)
	xlim2=c(56990716,57020327)

	for (pop in subpop) {

		sex=c('all','Homme','Femme')


		for (s in sex) {

			if(plot==1 || plot==2){
				s='all'
				ylim=c(0,0.15)
			}
			if(plot==3 || plot==4){
				s='Homme'
				ylim=c(0,0.15)
			}
			if(plot==5 || plot==6){
				s='Femme'
				ylim=c(0,0.25)
			}


			

			if(s=='all'){
				data_all<-read.table(paste(c("r2_",pop,".geno.ld"),collapse=''), header=TRUE, row.names=NULL,  sep="\t")
				
			}else{
				data_all<-read.table(paste(c("r2_",pop,".",s,".geno.ld"),collapse=''), header=TRUE, row.names=NULL,  sep="\t")
			}



			# Extract only our genes
			data_all=data_all[data_all$POS1>=xlim[1] & data_all$POS1<=xlim[2],]
			data_all=data_all[data_all$POS2>=xlim2[1] & data_all$POS2<=xlim2[2],]

			Couleur_6<-c("#56B4E9","#CC79A7","#D55E00","#009E73")

			# Get the 95th percentile
			percentile_95_All<-quantile(data_all[,5],0.95)

			# Extract our SNPs found in the Figure 3a
			data_all_7309=data_all[data_all$POS1==rs1967309,]
			data_CETP=data_all[data_all$POS2 %in% c(57007610,57008227,57008287,57014319),]


			data_all_7309$Color='black'
			data_all_7309[data_all_7309$POS2==57007610,'Color']=Couleur_6[1]
			data_all_7309[data_all_7309$POS2 %in% c(57008227),'Color']=Couleur_6[2]
			data_all_7309[data_all_7309$POS2 %in% c(57008287),'Color']=Couleur_6[3]
			data_all_7309[data_all_7309$POS2 == 57014319,'Color']=Couleur_6[4]

			
			data_CETP$Color='black'
			data_CETP[data_CETP$POS2==57007610,'Color']=Couleur_6[1]
			data_CETP[data_CETP$POS2 ==57008227,'Color']=Couleur_6[2]
			data_CETP[data_CETP$POS2 ==57008287,'Color']=Couleur_6[3]
			data_CETP[data_CETP$POS2 == 57014319,'Color']=Couleur_6[4]

			# CETP
			if(plot==1 || plot==3 || plot==5){
				plot(x=data_all_7309[,3]/diviseur, y=data_all_7309[,5], pch=16, xlab="Position in CETP (Mb)", ylab="", main=paste(c(pop,s),collapse='\n'), col=data_all_7309$Color,ylim=ylim,cex=2,yaxt='n', cex.axis=2)
				axis(2,ylim=range(data_all_7309[,5]),las=2,cex.axis=2)
				abline(h=percentile_95_All)
				legend('topright', c('rs158477','rs158480','rs158617','rs12447620'),col=Couleur_6,pch=16,bg="white",bty = "n",cex=2)
			}

			# ADCY9
			if(plot==2 || plot==4 || plot==6){
				plot(x=data_CETP[,2]/diviseur, y=data_CETP[,5], pch=16, xlab="Position in ADCY9 (Mb)", ylab="", main=paste(c(pop,unique(data_all[,4])),collapse='\n'), col=data_CETP$Color,ylim=ylim,cex=2,yaxt='n', cex.axis=2)
				axis(2,ylim=range(data_CETP[,5]),las=2,cex.axis=2)
				abline(h=percentile_95_All)
				abline(v=rs1967309/diviseur, lty=2,lwd=2)
				legend('topright', c('rs158477','rs158480','rs158617','rs12447620'),col=Couleur_6,pch=16,bg="white",bty = "n",cex=1.6)
				legend('topleft', c('rs1967309'),lty=2,lwd=2,bg="white",bty = "n",cex=1.5)
			}


			

		}
	}
	dev.off()
}