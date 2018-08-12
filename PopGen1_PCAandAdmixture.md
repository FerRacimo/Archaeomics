Exercises using PCA and Admixture
===============

# SSH

Let's log in to the servers (make sure you use the option -Y):

```
ssh -Y [YOUR KU ID]@ssh-snm-archaeo.science.ku.dk
ssh -Y archaeo-snm.science.domain
```

# Data and shortcuts

In this tutorial we will be using a combination of whole-genome and single-nucleotide polymorphism (SNP) capture datasets assembled by the Reich lab as part of a paleogenomic study of Europe and the Middle East (Lazaridis et al. 2016). The dataset contains genomic data from both ancient and present-day populations from around the world.

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
alias admixture='/science/groupdirs-nfs/SCIENCE-SNM-Archaeo/software/admixture_linux-1.3.0/admixture'
```

Now create a shortcut to your own data folder, and call that shortcut DATA. Thus, you will be easily able to dump all your intermediate data files. Type 'echo $DATA' to make sure you've created the shortcut successfully.

# LD pruning

We'll start with a set of files in plink format. Plink files usually come in sets of 3: a *fam* file, a *bim* file and a *bed* file. Our files will have the following names:

$HUMOR/AncientModern_reduced.bim

$HUMOR/AncientModern_reduced.fam

$HUMOR/AncientModern_reduced.bed

Take a look inside these three files. You can use the program *less*. What do you see inside? The bed file is compressed, and contains the genotypes of each individual in the *fam* file, ordered by the positions in the *bim* file.

When computing a PCA or performing an Admixture analysis, large datasets take a long time to analyze. However, a lot of SNPs actually have redundant information, as they may sit on the same haplotype and be in strong linkage disequilibrium (LD) with each other. We can "thin" our data to remove SNPs based on their LD correlation coefficients (r^2), using plink, keeping (almost) the same amount of information while significantly reducing the computational burden of our downstream algorithms. We can use the following commands to prune our data:

```
plink --bfile $HUMOR/AncientModern_reduced --indep-pairwise 50 10 0.1
plink --bfile $HUMOR/AncientModern_reduced --extract plink.prune.in --make-bed --out $DATA/AncientModern_reduced_pruned
```

The first command makes a list of SNPs that will be targeted for removal. These are SNPs with an r^2 value greater than 0.1 with any other SNP within a 50-SNP sliding window (with a 10-SNP overlap between windows). The second command performs the pruning.

Compare the number of SNPs in the $HUMOR/AncientModern_reduced.bim file and the $DATA/AncientModern_reduced_pruned.bim file. How many SNPs did we remove? How many did we keep?

# PCA

Now create a directory called 'PCA' and also create a shortcut for it called PCA.

We'll need to define a parameter file for convertf, which we'll call pca.par. Open your favorite text editor and write down the following lines. IMPORTANT: make sure your directory names are correctly written!

```
genotypename: [ YOUR DATA DIRECTORY NAME HERE ]/AncientModern_reduced_pruned.bed
snpname:      [ YOUR DATA DIRECTORY NAME HERE ]/AncientModern_reduced_pruned.bim
indivname:    [ YOUR DATA DIRECTORY NAME HERE ]/AncientModern_reduced_pruned.fam
evecoutname:  [ YOUR PCA DIRECTORY NAME HERE ]/eigenvectors.txt
evaloutname:  [ YOUR PCA DIRECTORY NAME HERE ]/eigenvalues.txt
numoutlieriter: 0
poplistname:  [YOUR PCA DIRECTORY NAME HERE]/pca_populations.txt
lsqproject: YES
```

Save this to a file called pca.par. The last two lines of this file server to ensure that the principal components will only be computed with the populations listed in the file "pca_populations.txt". As ancient DNA tends to be of lower quality and have more missing data than modern DNA, researchers traditionally calculate the principal components using modern genomes only, and then "project" the ancient genomes onto the modern PCA. Let's create the pca_popualtions.txt file:

 ```
 echo -e "Ju_hoan_North\nSardinian\nFrench\nItalian_North\nHan\nAmi\nYoruba\nMbuti\nPapuan\nOrcadian\nMayan\nKaritiana" > $PCA/pca_populations.txt
 ```

There's one last thing we need to do. For some weird reason, smartpca only recognizes population names if they're in the sixth column of the plink file, so we'll have to format the sixth column appropriately:

```
paste <(cat $DATA/AncientModern_reduced_pruned.fam | awk '{print $1,$2,$3,$4,$5,$1}') > temp.fam
mv temp.fam $DATA/AncientModern_reduced_pruned.fam
```

We can now run our PCA analysis using the smartpca program from eigensoft:

```
smartpca -p pca.par
```

To visualize the first 2 principal components from the PCA, we'll use an R script:

```
Rscript $SCRIPTS/plotPCA.R -a $PCA/eigenvalues.txt -e $PCA/eigenvectors.txt -c 1-2 -f $DATA/AncientModern_reduced_pruned.fam -o $PCA/PCA_World.pdf
xpdf $PCA/PCA_World.pdf
```

Which groups are separated along the first component of variation? Which groups are separated along the second component? Why do you think this is? How much of the total variance is captured by the first and second principal components?

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

We can visualize the admixture components of each individual using a barplot in R. For example, for K=5, we can type the following in R (make sure the directory name is correct):

```
R
K=5
Admdir <- [YOUR ADMIXTURE DIRECTORY HERE]
Datadir <- [YOUR DATA DIRECTORY HERE ] 
tbl <- read.table(paste(Admdir,"/AncientModern_reduced_pruned.",K,".Q",sep=""),header=FALSE)
inds <- read.table(paste(Datadir,"/AncientModern_reduced_pruned.fam",sep=""),header=FALSE) 
barplot(t(as.matrix(tbl)), col=rainbow(K), xlab="Individual #", ylab="Ancestry", border=NA,las=2,cex.names=0.3,names=inds[,1])
```

There are a lot of individuals here. We can also try plotting just a subset of them (say, from the 150th to the 250th):

```
indstoplot <- seq(150,250)
barplot(t(as.matrix(tbl))[,indstoplot], col=rainbow(K), xlab="Individual #", ylab="Ancestry", border=NA,las=2,cex.names=0.8,names=inds[indstoplot,1])
```

Do you notice any admixed individuals? How are Mayans modeled under this value of K?

Note that there are visualization programs better suited for exploring results from Admixture and other programs based on the Structure algorithms. One that is very useful and practical to use is Pong: https://github.com/ramachandran-lab/pong An R package that can also help with visualization of Admixture results is pophelper: https://github.com/royfrancis/pophelper

We can also look at the cross-validation errors:


```
grep -h CV log*.out
```

If you have experience with R, try to plot these values (otherwise work with a partner who does). Have they reached a minimum? If so, at which value of K? What does this mean?
