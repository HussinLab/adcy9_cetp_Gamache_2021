
import peer
import scipy as SP
import pandas as pd
import numpy as np
import os
import os.path
import sys


#Need file row=gene, col=ind
fil=sys.argv[1] # Individual tissue like this : 'Adipose-Subcutaneous'


file_to_analyse=('../Data/Expression/OnlyWhite/Normalized/'+fil+'.Gene.Filter6_20.OnlyWhite.normalizedLimmaVoom.txt')
file_output=('../Data/Expression/OnlyWhite/Normalized/'+fil+'.OnlyWhite.PEER.txt')


if os.path.isfile(file_to_analyse) :# and not os.path.isfile(file_output) :
	print(fil)

 	print("PEER: loading expression data ... ")
 	data = SP.loadtxt(file_to_analyse, delimiter="\t", dtype='string', comments='')
 	print("done.")



 	print("Modify data file")
 	# Get only expression
 	name_ind=data[0,1:]
 	name_gene=data[1:,0]
 	expr_file = data[1:,1:]
 	# Switch columns and rows
 	expr_file = np.ndarray.transpose(expr_file)
 	# Convert string to double
 	expr_file=[[float(y) for y in x] for x in expr_file]
	
 	print("done.")
 	print(name_gene)

	"""
		Let's begin!!
	"""

	#Need to adjust depending on the number of sample
 	nb_col=len(name_ind)
 	print(nb_col)
 	if (nb_col<16) :
	 	nb_hidden_factor=nb_col
	elif (nb_col<150) :
	 	nb_hidden_factor=15
	elif (nb_col<250):
		nb_hidden_factor=30
	elif (nb_col<350):
		nb_hidden_factor=45
	elif (nb_col>=350):
		nb_hidden_factor=60

 	print(nb_hidden_factor)



 	print("PEER: estimating hidden confounders ("+str(nb_hidden_factor)+")")

 	model = peer.PEER() # Create the model object
 	model.setPhenoMean(expr_file) #Set the observed data
 	model.setNk(nb_hidden_factor) # Infer K hidden confounders
	model.setNmax_iterations(100000) 
 	model.update() #Perform the inference

	print('Finish Parameters')	

 	factors=model.getX() #samples x PEER factors
 	precision = model.getAlpha() # PEER factors x 1
 	residuals=(model.getResiduals()) #samples x gene

 	print("Writing covariates file...")
 	inf_cov=["InferredCov" + str(i) for i in xrange(1,nb_hidden_factor+1)]
 	sam=["Sample_id"]
 	sam.extend(inf_cov)#Titre
 	sam=np.asmatrix(sam) #Convertis la liste en matrice pour fusion
 	df = pd.DataFrame(name_ind) #convertis liste des noms en matrice
 	frames=np.concatenate((df, factors), axis=1)# Fusionne les colonnes
 	temp=factors
 	frames=np.concatenate((sam, frames), axis=0)# Ajoute l'en-tete
 	np.savetxt(file_output, frames, fmt="%s", delimiter='\t', newline='\r\n')
 	print("done.")

 print('All data are ready for a simple linear regression with freely disponible covariates')
