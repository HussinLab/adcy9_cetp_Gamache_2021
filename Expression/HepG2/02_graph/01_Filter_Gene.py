import numpy as np
import scipy as SP

print('After script in 01_Filters')
print('Remove gene with low expresison')
geneIso='genes' # 

print('Read file 1')

# count file after removing ADCY9-205
data=np.loadtxt(("../01_Filters/ADCY9_Tardif.KO.GRCh38.E87.no205."+geneIso+".txt"), delimiter="\t", dtype='str')

ncol=data.shape[1]
nrow=data.shape[0]

keep=[0]
u=1
while u<nrow:
	i=1
	compte_higher=0
	while i<ncol:
		print(data[u,i])
		higher=(float(data[u,i])>=6) #Au moins 6 read
		compte_higher=compte_higher+higher
		i=i+1

	if compte_higher>=5:#In one echantillon
		keep.append(u)
	u=u+1
data1=data[keep,]
	

	
print('Write file')

np.savetxt(("RSEM_Hg38/"+geneIso+"ExpressionTable.Hg38.gene6_5.txt"),data1, fmt="%s", delimiter='\t', newline='\r\n')
	
