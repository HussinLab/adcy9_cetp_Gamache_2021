pheno=SEX
# file=SamplesInfo
file=Phenotypes

# Files from GTEx
n=$(grep -wn $pheno ../Data/Order$file.2.txt | cut -f 1 -d :)
echo $n
n2=1,$n
echo $n2

e=$(cut -f $n2 ../Data/FullFile/$file.txt)

echo $e >../Data/SpecificPhenotype/$pheno.txt

e1=$(cut -f 1 ../Data/FullFile/$file.txt)
e2=$(cut -f $n ../Data/FullFile/$file.txt)
echo $e1 >temp.txt
echo $e2 >>temp.txt

awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' temp.txt >../Data/SpecificPhenotype/$pheno.txt
head ../Data/SpecificPhenotype/$pheno.txt
