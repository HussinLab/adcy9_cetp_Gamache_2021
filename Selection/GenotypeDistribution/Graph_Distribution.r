#!/usr/bin/env Rscript
print('Figure 2a and 3c and suppl')

# The order in the graph
subpop2=c('AFR','YRI','LWK','GWD','MSL','ESN','ASW','ACB','EUR','CEU','TSI','FIN','GBR','IBS','EAS','CHB','JPT','CHS','CDX','KHV','SAS','GIH','PJL','BEB','STU','ITU','AMR','MXL','PUR','CLM','PEL','Native Americans','altaic','urralicYukaghir','chukchiKamchatkan','eskimoAleut','naDene','northernAmerind','centralAmerind','chibchanPaezan','equatorialTucanoan','andean')

# rs1967309 or rs158477
snp<-c("rs158477") 

if(snp=='rs158477'){
	name=c("GG","AG","AA")
	col=c('#56B4E9','#F0E442','#CC79A7')
	pos=57007610
}else if(snp=='rs1967309'){
	name=c("AA","AG","GG")
	col=c('#CC79A7','#F0E442','#56B4E9')
	pos=4065583
}

# Since NatAm are in a different files from 1000G
nat=FALSE


i<-1
all_data<-c()
# Create a table with all data
while(i<=length(subpop2)){
	if(subpop2[i]=='Native Americans'){
		nat=TRUE
	}
	if(nat){
		name<-paste(c("FN/hwe.",".hwe"), collapse=subpop2[i])
	}else{
		name<-paste(c("1000G/SNP_7309_Top_Geno_",".hwe"), collapse=subpop2[i])
	}
	if(!file.exists(name)){
		all_data<-c(all_data, 0,0,0)
	}else{
		data<-read.table(name, header=TRUE, sep="\t")
		print((data))
		if(sum(data$POS==pos)==0){
			all_data<-c(all_data, 0,0,0)
		}else{
			r1<-do.call(rbind, strsplit(as.character(data[data$POS==pos,3]), ('/')))
			data<-matrix(as.numeric(r1), ncol=3)
			tot<-data[1,1]+data[1,2]+data[1,3]
			all_data<-c(all_data, data[1,1], data[1,2], data[1,3])
		}
	}	
	i<-i+1
}


data_y<-matrix(all_data, nrow=3)

if(snp=='rs1967309'){
	pdf(paste(c(snp,'.pdf'),collapse=''),width=16, height=7)
}else{
	pdf(paste(c(snp,'.pdf'),collapse=''),width=16, height=7.5)
}


# Remove Bip pop names
data_x=c('','YRI','LWK','GWD','MSL','ESN','ASW','ACB','','CEU','TSI','FIN','GBR','IBS','','CHB','JPT','CHS','CDX','KHV','','GIH','PJL','BEB','STU','ITU','','MXL','PUR','CLM','PEL','','ALT','URY','CHK','ESA','NAD','NOA','CEA','CHP','EQT','AND')


# Remove pop with less than 40 samples
data_x=data_x[!(colSums(data_y)>0 & colSums(data_y)<40)]
data_y=data_y[,!(colSums(data_y)>0 & colSums(data_y)<40)]

# Add number to the names
for (i in seq(1, length(data_x), 1)) {
	if(sum(data_y[,i])!=0){
		s=sum(data_y[,i])
		data_x[i]=paste(c(data_x[i],' (',s,')'),collapse='')
	}
}

# Calcul proportion
data_y=prop.table(data_y,2)

par(mar=c(10, 5, 5, 2))

data_y=data_y[,-1]
data_x=data_x[-1]

barplot(data_y, main=snp, names.arg=data_x, col=col,xlab="", las=2, cex.names =1.9,cex=2, cex.main=2)

