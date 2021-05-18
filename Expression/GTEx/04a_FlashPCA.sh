#!/bin/sh


vcftools --gzvcf /project/6009524/shared/GTex/from_TerraBio/Genotypes/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files_GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz --remove-indels --min-alleles 2 --max-alleles 2 --max-missing 0.95 --recode --recode-INFO-all --out ../../Data/PCA/GTEx.v8.binalleles.maxMissing95

plink2 --vcf ../../Data/PCA/GTEx.v8.binalleles.maxMissing95.recode.vcf --indep-pairwise 1000 50 0.05
plink2 --vcf ../../Data/PCA/GTEx.v8.binalleles.maxMissing95.recode.vcf --extract plink2.prune.in --make-bed --out ../../Data/PCA/GTEx.v8.binalleles.maxMissing95.pruned


/project/6005588/shared/bin/flashpca_x86-64 --bfile ../../Data/PCA/GTEx.v8.binalleles.maxMissing95.pruned
