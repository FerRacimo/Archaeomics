Exercises using PCA and Admixture
===============

In this tutorial we will be using ...

We'll begin by downloading a gentoype data file from Lazaridis et al. 2016, containing SNP capture data from a variety of present-day and ancient populations:

...

Let's create a folder where we'll dump all our data files:

```
mkdir Data
```

# Data processing

Before we can start our analysis, we'll need to process the data in a certain way.

We only want to work with autosomal SNPs, so we'll first make a record of the SNPs that are located in chrX and chrY, so we can get rid of them later.

```
cat NearEastPublic/HumanOriginsPublic2068.snp | tr -s " " | awk 'BEGIN{OFS="\t"}{if ($2 == "23" || $2 == "24") print}' > Data/toremove.snp
```

Now, let's convert the genotype file from "packed eigenstrat" format to packed ped format. We'll need to define a parameter file for convertf, which we'll call geno2plink.par. Open your favorite text editor and write down the following lines:


```
genotypename:   NearEastPublic/HumanOriginsPublic2068.geno
snpname:        NearEastPublic/HumanOriginsPublic2068.snp
indivname:      NearEastPublic/HumanOriginsPublic2068.ind
outputformat:    PACKEDPED
genotypeoutname: Data/HumanOriginsPublic2068.bed
snpoutname:      Data/HumanOriginsPublic2068.bim
indivoutname:    Data/HumanOriginsPublic2068.fam
badsnpname:      Data/toremove.snp
```

Now, run convertf:

```
convertf -p geno2plink.par
```


# Fix *fam file
paste <(cat NearEastPublic/HumanOriginsPublic2068.ind | awk '{print $3}') <( cat Data/HumanOriginsPublic2068.fam | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$
5,$6}') > temp.fam
mv temp.fam Data/HumanOriginsPublic2068.fam
# Make list of families to extract
echo -e "Ju_hoan_North\nSardinian\nFrench\nItalian_North\tHan\nAmi\nYoruba\nMbuti\nPapuan\nOrcadian\nMayan\nKaritiana" > Data/groups_to_keep.txt
# Extract populations
plink --bfile Data/HumanOriginsPublic2068 --keep-fam Data/groups_to_keep.txt --make-bed --out Data/HumanOriginsPublic2068_reduced

# LD pruning


# PCA - takes time
mkdir PCA
smartpca -p pca.par
# Visualize PCA
Rscript scripts/PlotPCA.R PCA/eigenvectors.txt Data/HumanOriginsPublic2068_reduced.fam PCA/PCA_World.pdf


# Admixture analysis - takes time
mkdir Admixture
cd Admixture
mkdir K3
# Run Admixture with K=3
cd K3; admixture ../../Data/HumanOriginsPublic2068_reduced.bed 3; cd ..

# Visualize results with pong
echo -e "K3_run1\t3\tK3/HumanOriginsPublic2068_reduced.3.Q" > filemap.txt
