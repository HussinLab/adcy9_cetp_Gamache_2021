library(ggplot2)
library(RColorBrewer)
library("data.table")
print('Figure 1b')


freadf <- function(...) return(as.data.frame(fread(...)))

# Supplementa graph?
supp=FALSE

exon=read.table('Exon.ADCY9.txt',header=TRUE)

getRecombinaisonMap=function(data, gene){

	subpop=subpop
	data2=c()
	o=1

	for (i in subpop) {
		print(i)
		file=(paste(c('Chr16/',i,'-16-final.txt'),collapse=""))
		if(file.exists(file)){
			table=read.table(file, header=TRUE)
			print(head(table))
			data=data[data$gene==gene,]
			adcy9_line=(table[,1]>=range(data$x)[1] & table[,1]<=range(data$x)[2]) 
			datat=table[adcy9_line,]
			colnames(datat)[c(1,2)]=c('Position','Rate')
			datat$pop=i
			# Will change the color in the graph
			if(o<=7){
				datat$bigpop='AFR'
			}else if(o<=12){
				datat$bigpop='EUR'
			}else if(o<=17){
				datat$bigpop='EAS'
			}else if(o<=22){
				datat$bigpop='SAS'
			}else{
				datat$bigpop='AMR'
			}
			data2=rbind(data2,datat)
		}
		o=o+1
	}
	return(data2)
}

# Assign color to population
myColors=c('#E69F00','#009E73','#56B4E9','#D55E00','#CC79A7')
names(myColors)=c('AFR','EUR','EAS','SAS','AMR')
colScale=scale_colour_manual(name='bigpop', values=myColors)

plotGraph=function(data, NamePop, percentil=TRUE){
	d=data
	a_l=d$gene=='ADCY9'
	a=d[a_l,]
	c_l=d$gene=='CETP'
	c=d[c_l,]

	if(!supp){
		recom=getRecombinaisonMap(d, 'ADCY9')
	}
	
	a$bigpop=factor(a$bigpop,levels=c('AFR','EUR','EAS','SAS','AMR'))
	
	
	if(!supp){
		a=a[a$iHS>=2 | a$iHS<=-2,]
		recom$bigpop=factor(recom$bigpop,levels=c('AFR','EUR','EAS','SAS','AMR'))		
		recom$pop=factor(recom$pop,subpop)
		recom$Position=recom$Position/1000000
		a$iHS=abs(a$iHS)
	}
	
	a$pop=factor(a$pop,subpop)
	
	# Bp to Mb
	a$x=a$x/1000000
	xlim_a=xlim_a/1000000
	exon[,1]=exon[,1]/1000000
	exon[,2]=exon[,2]/1000000

	if(supp){
		g1=ggplot(a, aes(x=x,y=iHS, col=bigpop))
		g1=g1+geom_point(size=5)+ggtitle('')+xlab('')+ylab('')+xlim(xlim_a)+ylim(range(c(a$iHS,4)))
		g1=g1+geom_hline(yintercept = c(-2,2), color='grey')
		g1=g1+ colScale+ geom_segment(aes(x = 4065583/1000000, yend = -Inf, xend = 4065583/1000000, y = Inf), colour='grey', alpha=0.4)+annotate(x=4065583/1000000-0.003, y=3,size=10, geom='text',label='rs1967309', angle=90, alpha=0.4)

		g1=g1+geom_segment(aes(x=4045116/1000000,yend = -Inf, xend = 4045116/1000000, y = Inf), colour='black', alpha=0.3)
		g1=g1+geom_segment(aes(x=4077442/1000000,yend = -Inf, xend = 4077442/1000000, y = Inf), colour='black', alpha=0.3)
		
		g1=g1+facet_wrap(~pop, ncol=4)
		g1=g1+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position='none', strip.text=element_text(size=50), axis.text=element_text(size=50),panel.background=element_rect(fill='white',color='black'))
	}else{
		
		g1=ggplot(a, aes(x=x,y=iHS, col=bigpop), alpha=0.8)
		g1=g1+geom_point(size=3)+xlab('')+ylab('')+xlim(xlim_a)+ theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position=c(0.93,0.8351), legend.title=element_blank(), text=element_text(size=30), panel.border=element_rect(colour='black', linetype='solid'), axis.title=element_text(size=25),axis.text=element_text(color='black', size=25))
		
		g1=g1+geom_line(data=recom, linetype='twodash',alpha=0.8, aes(x=Position, y=(Rate/50)+2, color=bigpop))
		g1=g1+scale_y_continuous(sec.axis = sec_axis(~(.-2)*50, name = ""))
		g1=g1+geom_segment(x=min(exon$Begin), xend=max(exon$End),y=1.875, yend=1.875, colour='blue')
		g1=g1+geom_rect(xmax=min(exon$Begin),xmin=4012650/1000000,ymin=1.85,ymax=1.9, colour='blue', fill='blue')

		g1=g1+geom_rect(xmax=4166186/1000000,xmin=max(exon$End),ymin=1.85,ymax=1.9, colour='blue', fill='blue')
		g1=g1+geom_rect(data=exon, inherit.aes=FALSE,mapping=aes(xmin=Begin,xmax=End,ymin=1.8,ymax=1.95),  fill='blue')
		g1=g1+ colScale+ geom_segment(aes(x = 4065583/1000000, yend = 2.8, xend = 4065583/1000000, y = 3), colour='black',arrow = arrow(length = unit(0.03, "npc")))+annotate(x=4065583/1000000, y=3.5,size=10, geom='text',label='rs1967309', angle=90)

	}
	print(g1)


}






if(supp){
	png('iHS.JC.Sup2.png',width=3000, height=5000)
}else{
	png('iHS.JC.png',width=1000, height=500)
}


subpop=c('YRI','LWK','GWD','MSL','ESN','ASW','ACB','CEU','TSI','FIN','GBR','IBS','CHB','JPT','CHS','CDX','KHV','GIH','PJL','BEB','STU','ITU','MXL','PUR','CLM','PEL')

all_ADCY9=c()
all_CETP=c()

all_ADCY9=read.table('iHS.JC.txt', sep='\t',header=TRUE)

sup2_l=all_ADCY9[,2]>=2
sup2_adcy9=all_ADCY9[sup2_l,]


xlim_a=range(all_ADCY9$x)
plotGraph(all_ADCY9,'all')



