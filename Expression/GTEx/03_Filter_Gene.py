import numpy as np
import scipy as SP
import os.path
from os import path

geneIso='Isoforms' # Genes 	Isoforms

print('Read file')
# File with all tissues
tissu=np.loadtxt('AllTissues.txt',  delimiter="\t", dtype='str')

print(tissu)

for x in tissu:
	print(x)
	print(path.isfile("../Data/Expression/Filtrer/"+geneIso+"/"+x+"."+geneIso+".Filter6_20.txt"))
	if path.isfile("../Data/Expression/Filtrer/"+geneIso+"/"+x+"."+geneIso+".Filter6_20.txt"):
		continue

	print('No pass')
	data=np.loadtxt(("../Data/Expression/Brut/"+geneIso+"/"+x+"."+geneIso+".gct"),  delimiter="\t", dtype='str', skiprows=2)

	data=np.delete(data,1,1)


	print('Count')

	ncol=data.shape[1]
	nrow=data.shape[0]

	# Filter
	keep=[0]
	u=1
	while u<nrow:
		i=1
		compte_higher=0
		while i<ncol:
			higher=(float(data[u,i])>=6) #Au moins 1 read
			compte_higher=compte_higher+higher
			i=i+1

		passing=float(compte_higher)/float(ncol-1)
		
		if passing>=0.2:
			keep.append(u)
		print(u)
		u=u+1
	data1=data[keep,]



	print('Write file')
	np.savetxt(("../Data/Expression/Filtrer/"+geneIso+"/"+x+"."+geneIso+".Filter6_20.txt"),data1, fmt="%s", delimiter='\t', newline='\r\n')

