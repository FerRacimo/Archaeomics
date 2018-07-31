Exercises using PCA and Admixture
===============

In this tutorial we will be using a combination of whole-genome and SNP capture datasets assembled by the Reich lab as part of a paleogenomic study of Europe and the Middle East (Lazaridis et al. 2016). The dataset contains genomic data from both ancient and present-day populations from around the world.

The data files can be found here:

```
/science/groupdirs-nfs/SCIENCE-SNM-Archaeo/Day3_data/HumanOriginsData
```

Let's create a shortcut to this folder, so that we can easily reference it when working from our own personal folders. We'll also create an alias for certain programs that we will use  below.

```
HUMOR=/science/groupdirs-nfs/SCIENCE-SNM-Archaeo/class_aDNA_2018/Day3_data/HumanOriginsData
SCRIPTS=/science/groupdirs-nfs/SCIENCE-SNM-Archaeo/class_aDNA_2018/scripts
alias plink='/science/groupdirs-nfs/SCIENCE-SNM-Archaeo/software/plink1.9/plink'
alias smartpca='/science/groupdirs-nfs/SCIENCE-SNM-Archaeo/software/EIG/bin/smartpca'
alias convertf='/science/groupdirs-nfs/SCIENCE-SNM-Archaeo/software/EIG/bin/convertf'
alias mergeit='/science/groupdirs-nfs/SCIENCE-SNM-Archaeo/software/EIG/bin/mergeit'
alias adimxture='/science/groupdirs-nfs/SCIENCE-SNM-Archaeo/software/admixture_linux-1.3.0/admixture'
```


Now create a shortcut to your own data folder, and call that shortcut DATA. Thus, you will be easily able to dump all your intermediate data files. Type 'echo $DATA' to make sure you've created the shortcut successfully.


# Data processing

Before we can start our analysis, we'll need to clean our data for downstream analyses. We only want to work with autosomal SNPs, so we'll first make a record of the SNPs that are located in chrX and chrY, so we can get rid of them later.

```
cat $HUMOR/AncientModern.snp | tr -s " " | awk 'BEGIN{OFS="\t"}{if ($2 == "23" || $2 == "24") print}' > $DATA/toremove.snp
```

Now, let's convert the genotype file from "packed eigenstrat" format to "packed ped format. We'll need to define a parameter file for convertf, which we'll call geno2plink.par. Open your favorite text editor and write down the following lines

IMPORTANT: you need to replace the DATA line for whichever Data folder you're working in!


```
HUMOR: /science/groupdirs-nfs/SCIENCE-SNM-Archaeo/class_aDNA_2018/Day3_data/HumanOriginsData 
DATA: [ YOUR DATA FOLDER HERE ]
genotypename:   HUMOR/AncientModern.geno
snpname:        HUMOR/AncientModern.snp
indivname:      HUMOR/AncientModern.ind
outputformat:    PACKEDPED
genotypeoutname: DATA/AncientModern.bed
snpoutname:      DATA/AncientModern.bim
indivoutname:    DATA/AncientModern.fam
badsnpname:      DATA/toremove.snp
```

Save the file and name it 'geno2plink.par'. Now, run convertf:

```
convertf -p geno2plink.par
```

We now need to fix the *fam file to remove unwanted spaces:

```
paste <(cat $HUMOR/AncientModern.ind | awk '{print $3}') <( cat $DATA/AncientModern.fam | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5,$6}') > temp.fam
mv temp.fam $DATA/AncientModern.fam
```

Now, we'll make a list of populations ("plink families") to focus on for donwstream analyses, and extract them from the plink file. We'll extract a subset of modern and ancient populations:

```
echo -e "Ju_hoan_North\nSardinian\nFrench\nItalian_North\tHan\nAmi\nYoruba\nMbuti\nPapuan\nOrcadian\nMayan\nKaritiana\nEurope_LNBA\nSteppe_EMBA" > $DATA/groups_to_keep.txt
plink --bfile $DATA/AncientModern --keep-fam $DATA/groups_to_keep.txt --make-bed --out $DATA/AncientModern_reduced
```

# LD pruning

When computing a PCA or performing an Admixture analysis, large datasets take a long time to analyze. However, a lot of SNPs actually have redundant information, as they may sit on the same haplotype and be in strong linkage disequilibrium with each other. We can "thin" our data to remove SNPs based on their LD correlation coefficients, using plink, keeping (almost) the same amount of SNP data while significantly reducing the computational burden of our downstream algorithms. We can use the following commands to prune our data:

```
plink --bfile $DATA/AncientModern_reduced --indep-pairwise 50 10 0.1
plink --bfile $DATA/AncientModern_reduced --extract plink.prune.in --make-bed --out $DATA/AncientModern_reduced_pruned

```

The first command makes a list of SNPs that will be targeted for removal. These are SNPs with an r^2 value greater than 0.1 with any other SNP within a 50-SNP sliding window (with a 10-SNP overlap between windows). The second command performs the pruning.

Compare the number of SNPs in the $DATA/AncientModern_reduced.bim file and the $DATA/AncientModern_reduced_pruned.bim file. How many SNPs did we remove? How many did we keep?

# PCA

Now create a directory called 'PCA' and also create a shortcut for it called PCA.

Let's create a parameter file for PCA (called pca.par), using your favorite text editor. IMPORTANT: make sure your directory names are correctly written!

```
genotypename: [ YOUR DATA DIRECTORY NAME HERE ]/AncientModern_reduced_pruned.bed
snpname:      [ YOUR DATA DIRECTORY NAME HERE ]/AncientModern_reduced_pruned.bim
indivname:    [ YOUR DATA DIRECTORY NAME HERE ]/AncientModern_reduced_pruned.fam
evecoutname:  [ YOUR PCA DIRECTORY NAME HERE ]/eigenvectors.txt
evaloutname:  [ YOUR PCA DIRECTORY NAME HERE ]/eigenvalues.txt
numoutlieriter: 0
```

Save this to a file called pca.par. We can now run our PCA analysis using the smartpca program from eigensoft:

```
smartpca -p pca.par
```

To visualize the first 2 principal components from the PCA, we'll use an R script:

```
Rscript $SCRIPTS/plotPCA.R -a $PCA/eigenvalues.txt -e $PCA/eigenvectors.txt -c 1-2 -f $DATA/AncientModern_reduced_pruned.fam -o $PCA/PCA_World.pdf
xpdf $PCA/PCA_World.pdf
```

Which groups are separated along the first component of variation? Which groups are separated along the second component? Why do you think this is?

Look at the file containing the eigenvalues. The value of a particular eigenvalue divided by the total sum of all eigenvalues is the proportion of variance explained by its corresponding eigenvector (principal component). Load the eigenvalue table into R, and calculate how much variance in the data is captured by the first, second, third and fourth principal components.

# Admixture analysis

Let's create a directory to dump our Admixture ouput files:

```
mkdir Admixture
```

Now, we'll run the Admixture program with K=1, K=2, K=3 and K=4 and K=5 ancestral components. Note that this will take some time to run as the algorithm performs a series of EM steps until it converges. In the mean time, you may want to review the lectures from today (or check facebook, your call). We'll use the --cv flag so that the program also performs 5-fold cross-validation in each run.

```
cd Admixture
for K in {1..5}; do
admixture --cv $DATA/AncientModern_reduced_pruned.bed $K  | tee log_${K}.out
done
cd ..
```

Look at the output files in the Admixture folder. Based on today's lecture, what do you think each of these represent?

We can visualize the admixture components of each individual using a barplot in R. For example, for K=3, we can type the following in R (make sure the directory names are correct):

```
R
K=3
tbl <- read.table(paste("Admixture/AncientModern_reduced_pruned.",K,".Q",sep=""),header=FALSE)
inds <- read.table("Data/AncientModern_reduced_pruned.fam",header=FALSE) 
barplot(t(as.matrix(tbl)), col=rainbow(K), xlab="Individual #", ylab="Ancestry", border=NA,las=2,cex.names=0.3,names=inds[,1])
```


Note that there are visualization programs better suited for exploring results from Admixture and other programs based on the Structure algorithms. One that is very useful and practical to use is Pong: https://github.com/ramachandran-lab/pong

We can also look at the cross-validation errors:


```
grep -h CV log*.out
```

If you have experience with R, try to plot these values (otherwise work with a partner who does). Have they reached a minimum? If so, at which value of K? What does this mean?
