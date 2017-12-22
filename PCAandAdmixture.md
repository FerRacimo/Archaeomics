Exercises using PCA and Admixture
===============

In this tutorial we will be using a SNP capture dataset produced by the Reich lab as part of a paleogenomic study of Europe and the Middle East (Lazaridis et al. 2016). The dataset contains genomes from both ancient and present-day populations.

We'll begin by downloading the genotype data file:

```
wget https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/NearEastPublic.tar.gz
tar xvzf NearEastPublic.tar.gz
```

Let's create a folder where we'll dump all our intermediate data files:

```
mkdir Data
```


# Data processing

Before we can start our analysis, we'll need to clean our data for downstream analyses. We only want to work with autosomal SNPs, so we'll first make a record of the SNPs that are located in chrX and chrY, so we can get rid of them later.

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

We now need to fix the *fam file to remove unwanted spaces:

```
paste <(cat NearEastPublic/HumanOriginsPublic2068.ind | awk '{print $3}') <( cat Data/HumanOriginsPublic2068.fam | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5,$6}') > temp.fam
mv temp.fam Data/HumanOriginsPublic2068.fam
```

Now, we'll make a list of populations ("plink families") to focus on for donwstream analyses, and extract them from the plink file:

```
echo -e "Ju_hoan_North\nSardinian\nFrench\nItalian_North\tHan\nAmi\nYoruba\nMbuti\nPapuan\nOrcadian\nMayan\nKaritiana" > Data/groups_to_keep.txt
plink --bfile Data/HumanOriginsPublic2068 --keep-fam Data/groups_to_keep.txt --make-bed --out Data/HumanOriginsPublic2068_reduced
```

# LD pruning

When computing a PCA or performing an Admixture analysis, large datasets take a long time to analyze. However, a lot of SNPs actually have redundant information, as they may sit on the same haplotype and be in strong linkage disequilibrium with each other. We can "thin" our data to remove SNPs based on their LD correlation coefficients, using plink, keeping (almost) the same amount of SNP data while significantly reducing the computational burden of our downstream algorithms. We can use the following commands to prune our data:

```
plink --bfile Data/HumanOriginsPublic2068_reduced --indep-pairwise 50 10 0.1
plink --bfile Data/HumanOriginsPublic2068_reduced --extract plink.prune.in --make-bed --out Data/HumanOriginsPublic2068_reduced_pruned

```

The first command makes a list of SNPs that will be targeted for removal. These are SNPs with an r^2 value greater than 0.1 with any other SNP within a 50-SNP sliding window (with a 10-SNP overlap between windows). The second command performs the pruning.

Compare the number of SNPs in the Data/HumanOriginsPublic2068_reduced.bim file and the Data/HumanOriginsPublic2068_reduced_pruned.bim file. How many SNPs did we remove? How many did we keep?

# PCA

Let's create a directory for our PCA results:
```
mkdir PCA
```

Now, create a parameter file for PCA (called pca.par), using your favorite text editor:

```
genotypename: Data/HumanOriginsPublic2068_reduced_pruned.bed
snpname:      Data/HumanOriginsPublic2068_reduced_pruned.bim
indivname:    Data/HumanOriginsPublic2068_reduced_pruned.fam
evecoutname:  PCA/eigenvectors.txt
evaloutname:  PCA/eigenvalues.txt
numoutlieriter: 0
```

We can now run our PCA analysis using the smartpca program from eigensoft:

```
smartpca -p pca.par
```

To visualize the first 2 principal components from the PCA, we'll use an R script:

```
Rscript scripts/PlotPCA.R PCA/eigenvectors.txt Data/HumanOriginsPublic2068_reduced_pruned.fam PCA/PCA_World.pdf
xpdf PCA/PCA_World.pdf
```

Which groups are separated along the first component of variation? How much variance in the data is captured by this component? Which groups are separated along the second component? Why do you think this is?

# Admixture analysis

```
mkdir Admixture
cd Admixture
mkdir K3
```

Run Admixture with K=3

```
cd K3; admixture ../../Data/HumanOriginsPublic2068_reduced_pruned.bed 3; cd ..

```

Visualize results with pong

```
echo -e "K3_run1\t3\tK3/HumanOriginsPublic2068_reduced_pruned.3.Q" > filemap.txt

```
