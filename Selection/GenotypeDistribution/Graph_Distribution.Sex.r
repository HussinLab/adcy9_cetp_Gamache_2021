#!/usr/bin/env Rscript
print('Figure 2a and 3c, Figure 4-figure supplement 1')


# The order in the graph
subpop2=c('AFR','YRI','LWK','GWD','MSL','ESN','ASW','ACB','EUR','CEU','TSI','FIN','GBR','IBS','EAS','CHB','JPT','CHS','CDX','KHV','SAS','GIH','PJL','BEB','STU','ITU','AMR','MXL','PUR','CLM','PEL','Native Americans','altaic','urralicYukaghir','chukchiKamchatkan','eskimoAleut','naDene','northernAmerind','centralAmerind','chibchanPaezan','equatorialTucanoan','andean')

snp<-c("rs1967309") 
# Option : rs1967309 or rs158477

# If stratified by sex?
bySex=FALSE

if(bySex){
	subpop2=c(subpop2,'LIMAA2','LIMAA')
}

sex='MALE'

if(snp=='rs158477'){
	name=c("GG","AG","AA")
	col=c('#56B4E9','#F0E442','#CC79A7')
	pos=57007610
}else if(snp=='rs1967309'){
	name=c("AA","AG","GG")
	col=c('#CC79A7','#F0E442','#56B4E9')
	pos=4065583
}else{
	print('snp unknown - error')
	print()
}

# Since NatAm/LIMAA are in a different files from 1000G
nat=FALSE #NAGD?
lim=FALSE #LIMAA?


i<-1
ALL<-c()

# Get all info in one table
while(i<=length(subpop2)){
	if(subpop2[i]=='Native Americans'){
		nat=TRUE

	}
	if(subpop2[i]=='LIMAA'){
		nat=FALSE
		lim=TRUE
	}

	if(bySex){
		if(nat){
			if(bySex){
				if(sex=='FEMALE' || sex=='Females'){
					sex='Females'
				}else{
					sex='Males'
				}
			}
			Name<-paste(c("FN/hwe.",subpop2[i],'.',sex,".hwe"), collapse='')
		}else if(lim){
			if(bySex){
				if(sex=='FEMALE' || sex=='Females'){
					sex='Female'
				}else{
					sex='Male'
				}
			}
			Name<-paste(c(subpop2[i],'/SNPs2.',sex,".hwe"), collapse='')
		}else{
			Name<-paste(c("1000G/ADCY9_CETP.",subpop2[i],'.',sex,".rs1967309.rs158477.hwe"), collapse='')
		
		}
	}else{
		if(nat){
			Name<-paste(c("FN/hwe.",".hwe"), collapse=subpop2[i])
		}else if(lim){
			Name<-paste(c(subpop2[i],'/SNPs2.All.hwe'), collapse='')
		}else{
			Name<-paste(c("1000G/SNP_7309_Top_Geno_",".hwe"), collapse=subpop2[i])
		}
	}
	

	if(!file.exists(Name)){
		# If it is for a big pop, we put zero to be visual friendly 
		print(Name)
		ALL<-c(ALL, 0,0,0)
	}else{
		Data<-read.table(Name, header=TRUE, sep="\t")
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


data_x<-subpop2
data_y<-matrix(ALL, nrow=3)


if(!bySex){
	pdf(paste(c(snp,'.pdf'),collapse=''),width=16, height=7)
	
}else{
	
	pdf(paste(c(snp,'.',sex,'.pdf'),collapse=''),width=16, height=7)
}


# Remove Bip pop names
data_x=c('','YRI','LWK','GWD','MSL','ESN','ASW','ACB','','CEU','TSI','FIN','GBR','IBS','','CHB','JPT','CHS','CDX','KHV','','GIH','PJL','BEB','STU','ITU','','MXL','PUR','CLM','PEL','','ALT','URY','CHK','ESA','NAD','NOA','CEA','CHP','EQT','AND')

if(bySex){
	data_x=c(data_x,'','LIMAA')
}


# Remove pop with too few individual
if(sex=='Males'){
	data_x=data_x[!(colSums(data_y)>0 & colSums(data_y)<13)]
	data_y=data_y[,!(colSums(data_y)>0 & colSums(data_y)<13)]
	par(mar=c(11, 5, 5, 2))
}else if(sex=='Females'){
	data_x=data_x[!(colSums(data_y)>0 & colSums(data_y)<13)]
	data_y=data_y[,!(colSums(data_y)>0 & colSums(data_y)<13)]
	# To get the same pop than Male
	data_y=data_y[,data_x!='ESA']
	data_x=data_x[data_x!='ESA']
	par(mar=c(11, 5, 5, 2))
}else{
	data_x=data_x[!(colSums(data_y)>0 & colSums(data_y)<40)]
	data_y=data_y[,!(colSums(data_y)>0 & colSums(data_y)<40)]
	par(mar=c(10, 5, 5, 2))
}


# Add number to the names
for (i in seq(1, length(data_x), 1)) {
	if(sum(data_y[,i])!=0){
		s=sum(data_y[,i])
		data_x[i]=paste(c(data_x[i],' (',s,')'),collapse='')
	}
}

# Calcul proportion
data_y=prop.table(data_y,2)



data_y=data_y[,-1]
data_x=data_x[-1]

barplot(data_y, main=snp, names.arg=data_x, col=col,xlab="", las=2, cex.names =1.9,cex=2, cex.main=2)

