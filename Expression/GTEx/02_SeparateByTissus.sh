
# Separe the count file by tissues
while read t
do
	echo $t
	a=$(grep -nf ../Data/Covariables/ID.PassFilter10Millions.MaxDuplicated.$t.txt ../Data/OrdreIDExpressionFiles.txt | cut -f 1 -d : | tr '\n' ',')

	b='1,2,'$a
	b=${b::-1}

	# Genes
	c=$(cut -f $b /project/6009524/shared/GTex/from_TerraBio/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct >../Data/Expression/Brut/Genes/$t.Gene.gct)
done<AllTissues.txt



