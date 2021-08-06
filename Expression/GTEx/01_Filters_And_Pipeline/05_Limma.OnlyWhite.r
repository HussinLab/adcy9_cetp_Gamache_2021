#!/usr/bin/env Rscript
library("data.table")
library("plyr")
library("limma")
library("edgeR")
library("reshape2")
library("ggrepel")
library("ggfortify")
library('ggplot2')

freadf <- function(...) return(as.data.frame(fread(...)))

# List of all tissues
organlist=read.table('AllTissu.txt')


r=1
while(r<=nrow(organlist)){
	organ=as.character(organlist[r,1])
	print(organ)
	print("Read reads")


	if(!file.exists(paste(c('../Data/Expression/OnlyWhite/Normalized/',organ,'.Gene.Filter6_20.OnlyWhite.normalizedLimmaVoom.txt'),collapse=""))){

		print('Will need to extract individual obtained from 04b from file from 03')
		mydata<-freadf(paste(c('../Data/Expression/OnlyWhite/',organ,'.Gene.Filter6_20.OnlyWhite.txt'),collapse=""))

		if(ncol(mydata)>10){

			name<-colnames(mydata)
			rownames(mydata)<-mydata[,1]#Met le numéro du transcript comme nom de ligne
			gene<-as.data.frame(mydata[,1])
			mydata<-mydata[,-c(1)]


			library("dplyr")
			print("Normalisation avec voom (JC)")


			md<-DGEList(mydata) 


			# Convert TPM values to log cpm
			lcpm <- cpm(md, log=TRUE)
			dge<-md



			melt_lcpm <- melt(lcpm)


			median_melt_lcpm <- melt_lcpm %>%
				group_by(Var1) %>%
				summarize(value = median(value))


			## Thresholding: filter out lowly expressed genes that could throw off the analyses
			# Keep genes for which the mean lcpm value is above -5
			average_lcpm <- aveLogCPM(dge)
			keep_lcpm<-((average_lcpm>-5))

			dge<-dge[keep_lcpm,]
			gene<-gene[keep_lcpm,]

			# Show the effect of normalization
			raw_dge <- dge
			norm_dge <- calcNormFactors(dge, method = "TMM")

			raw_lcpm <- cpm(raw_dge, log = TRUE)
			norm_lcpm <- cpm(norm_dge, log = TRUE)


			#Normalize with voom
			voomed<-voom(norm_dge, plot=T)
			voomed_reads<-as.data.frame(voomed$E)
			name_row<-as.matrix(rownames(voomed_reads))
			colnames(name_row)<-c("Name")
			voomed_reads<-cbind(name_row,voomed_reads)
			norm_lcpm<-(cbind(name_row,norm_lcpm))


			fwrite(voomed_reads,file=paste(c('../Data/Expression/OnlyWhite/Normalized/',organ,'.Gene.Filter6_20.OnlyWhite.normalizedLimmaVoom.txt'),collapse=""), quote=FALSE, sep="\t", col.names=TRUE)
		}
	}
	r=r+1
}



	
