library('reshape2')

# Info extracted from information file from GTEx, two columns : SAMPID et SMMPPD/SMTSD

# Coverage by ID
info=read.table('../../../Phenotype/Phenotype_GTEx/Data/SpecificPhenotype/SMMPPD.txt', fill=TRUE, header=TRUE)

# Association ID-tissue
tissu=read.table('../../../Phenotype/Phenotype_GTEx/Data/SpecificPhenotype/SMTSD.txt', fill=TRUE, header=TRUE)

all=merge(info,tissu,by='SAMPID')

a=colsplit(string=tissu[,1],pattern='-',names=c('A','B','C','D','E'))
all$ID=a$B

# Filter on coverage and associated to a tissue
higher10millions=all$SMMPPD>10000000 & !is.na(all$SMMPPD)
pass=all[higher10millions,]

allTissu=unique(pass$SMTSD)

write.table(allTissu,'AllTissues.txt',row.names=FALSE,col.names=FALSE)
keep=c()

# For individual who are duplicate, keep the one with the highest coverage
for (i in allTissu) {
	print(i)
	thisTissu=pass$SMTSD==i
	oneID=unique(pass[thisTissu,'ID'])

	subKeep=c()
	for (j in oneID) {
		thisID=pass$ID==j
		duplicatedInThisTissu=thisTissu & thisID

		if(sum(duplicatedInThisTissu)>1){
			
			m=(which.max(pass[duplicatedInThisTissu,2]))
			p=pass[duplicatedInThisTissu,]
			m=p[m,1]
			keep=c(keep,as.character(m))
			subKeep=c(subKeep,as.character(m))

		}else{
			keep=c(keep,as.character(pass[duplicatedInThisTissu,1]))
			subKeep=c(subKeep,as.character(pass[duplicatedInThisTissu,1]))

		}

	}
	passDup=pass$SAMPID %in% subKeep
	pass2=pass[passDup,]

	idOnly=as.data.frame(as.character(pass2[,1]))
	write.table(idOnly,paste(c('../Data/Covariables/ID.PassFilter10Millions.MaxDuplicated.',i,'.txt'),collapse=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
}

passDup=pass$SAMPID %in% keep
pass2=pass[passDup,]

idOnly=as.data.frame(as.character(pass2[,1]))

# Write ID who pass
write.table(idOnly,'../Data/Covariables/ID.PassFilter10Millions.MaxDuplicated.txt',quote=FALSE,row.names=FALSE,col.names=FALSE)
