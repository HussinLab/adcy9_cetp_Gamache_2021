#!/usr/bin/env Rscript
print('Figure 1a and 2c')

# The order in the graph
subpop2=c('AFR','YRI','LWK','GWD','MSL','ESN','ASW','ACB','EUR','CEU','TSI','FIN','GBR','IBS','EAS','CHB','JPT','CHS','CDX','KHV','SAS','GIH','PJL','BEB','STU','ITU','AMR','MXL','PUR','CLM','PEL','Native Americans','altaic','urralicYukaghir','chukchiKamchatkan','eskimoAleut','naDene','northernAmerind','centralAmerind','chibchanPaezan','equatorialTucanoan','andean')

# rs1967309 or rs158477
SNP<-c("rs158477") 

if(SNP=='rs158477'){
	name=c("GG","AG","AA")
	col=c('#56B4E9','#F0E442','#CC79A7')
	pos=57007610
}else if(SNP=='rs1967309'){
	name=c("AA","AG","GG")
	col=c('#CC79A7','#F0E442','#56B4E9')
	pos=4065583
}

# Since NatAm are in a different files from 1000G
nat=FALSE
v<-1


i<-1
ALL<-c()
# Create a table with all data
while(i<=length(subpop2)){
	if(subpop2[i]=='Native Americans'){
		nat=TRUE
	}
	if(nat){
		Name<-paste(c("FN/hwe.",".hwe"), collapse=subpop2[i])
	}else{
		Name<-paste(c("1000G/SNP_7309_Top_Geno_",".hwe"), collapse=subpop2[i])
	}
	if(!file.exists(Name)){
		ALL<-c(ALL, 0,0,0)
	}else{
		Data<-read.table(Name, header=TRUE, sep="\t")
		print((Data))
		if(sum(Data$POS==pos)==0){
			ALL<-c(ALL, 0,0,0)
		}else{
			R1<-do.call(rbind, strsplit(as.character(Data[Data$POS==pos,3]), ('/')))
			Data<-matrix(as.numeric(R1), ncol=3)
			tot<-Data[1,1]+Data[1,2]+Data[1,3]
			ALL<-c(ALL, Data[1,1], Data[1,2], Data[1,3])
		}
	}	
	i<-i+1
}

Donnex<-subpop2
Donne1y<-matrix(ALL, nrow=3)

png(paste(c(SNP,'.png'),collapse=''),width=1000, height=500)

# Remove Bip pop names
Donnex=c('','YRI','LWK','GWD','MSL','ESN','ASW','ACB','','CEU','TSI','FIN','GBR','IBS','','CHB','JPT','CHS','CDX','KHV','','GIH','PJL','BEB','STU','ITU','','MXL','PUR','CLM','PEL','','ALT','URY','CHK','ESA','NAD','NOA','CEA','CHP','EQT','AND')


density1 <- seq(45,50,length.out=1)

# Remove pop with less than 40 samples
Donnex=Donnex[!(colSums(Donne1y)>0 & colSums(Donne1y)<40)]
Donne1y=Donne1y[,!(colSums(Donne1y)>0 & colSums(Donne1y)<40)]

# Add number to the names
for (i in seq(1, length(Donnex), 1)) {
	if(sum(Donne1y[,i])!=0){
		s=sum(Donne1y[,i])
		Donnex[i]=paste(c(Donnex[i],' (',s,')'),collapse='')
	}
}

# Calcul proportion
Donne1y=prop.table(Donne1y,2)

par(mar=c(10, 5, 5, 2))

Donne1y=Donne1y[,-1]
Donnex=Donnex[-1]

barplot(Donne1y, main=SNP[v], names.arg=Donnex, col=col,xlab="", las=2, cex.names =1.9,cex=2, cex.main=2)
legend(par('usr')[1]+1,par('usr')[4]+0.2, xpd=TRUE, name, fill= col, ncol=3, bty = "n", cex=2)
