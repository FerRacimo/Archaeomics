
# Data processing

Before we start, let's define some shortcuts again:

```
HUMOR=/science/groupdirs-nfs/SCIENCE-SNM-Archaeo/class_aDNA_2018/Day3_data/HumanOriginsData
SCRIPTS=/science/groupdirs-nfs/SCIENCE-SNM-Archaeo/class_aDNA_2018/scripts
alias plink='/science/groupdirs-nfs/SCIENCE-SNM-Archaeo/software/plink1.9/plink'
alias convertf='/science/groupdirs-nfs/SCIENCE-SNM-Archaeo/software/EIG/bin/convertf'
alias qp3Pop='/science/groupdirs-nfs/SCIENCE-SNM-Archaeo/software/AdmixTools/bin/qp3Pop'
alias qpDstat='/science/groupdirs-nfs/SCIENCE-SNM-Archaeo/software/AdmixTools/bin/qpDstat'
```

You will also need to re-define the alias to your particular $DATA folder.

Now, let's convert our pruned plink files to eigenstrat format. Write a new parameter file (let's call it plink2geno.par), with the following information:

```
genotypename: [YOUR DATA DIRECTORY HERE]/AncientModern_reduced_pruned.bed
snpname: [YOUR DATA DIRECTORY HERE]/AncientModern_reduced_pruned.bim
indivname: [YOUR DATA DIRECTORY HERE]/AncientModern_reduced_pruned.fam
outputformat: EIGENSTRAT
genotypeoutname: [YOUR DATA DIRECTORY HERE]/AncientModern_reduced_pruned.geno
snpoutname: [YOUR DATA DIRECTORY HERE]/AncientModern_reduced_pruned.snp
indivoutname: [YOUR DATA DIRECTORY HERE]/AncientModern_reduced_pruned.ind
```

Run convertf:

```
convertf -p plink2geno.par
```

Before we proceed, we need to clean up our ind file:

```
paste <(cut -d ":" -f 2 $DATA/AncientModern_reduced_pruned.ind | cut -d " "  -f 1) <(cut -d ":" -f 2 $DATA/AncientModern_reduced_pruned.ind | cut -d " " -f 2) <(cut -d ":" -f 1  $DATA/AncientModern_reduced_pruned.ind | sed 's/ //g') > $DATA/temp.txt
mv $DATA/temp.txt $DATA/AncientModern_reduced_pruned.ind
```

Take a look inside the new *ind*, *geno* and *snp* files (you can use the program *less*). What do you see inside? How are they different from the plink format files?

# Outgroup F3 statistics

We'll begin by computing outgroup F3 statistics to determine which are the populations that share the most drift with French. We'll use Ju_hoan_North - a population from Namibia - as the outgroup. A requirement when computing F3 statistics in AdmixTools is to have a list of population triplets. Copy the following list of triplets to a file and name it OutgroupF3_French.txt. Note that the first column corresponds to our population of interest (French), the second column iterates over a series of populations we'll compare the French to, and the third column is a fixed outgroup, in this case Ju_hoan_North.

```
French  Papuan  Ju_hoan_North
French  Sardinian Ju_hoan_North
French  Orcadian  Ju_hoan_North
French  Mayan Ju_hoan_North
French  Yoruba  Ju_hoan_North
French  Mbuti Ju_hoan_North
French  Karitiana Ju_hoan_North
French  Italian_North Ju_hoan_North
French  Ami Ju_hoan_North
```

We also need to create a parameter file (OutgroupF3.par) to specify the location of our input files:

```
genotypename:   [YOUR DATA DIRECTORY HERE]/AncientModern_reduced_pruned.geno
snpname:	[YOUR DATA DIRECTORY HERE]/AncientModern_reduced_pruned.snp
indivname:	[YOUR DATA DIRECTORY HERE]/AncientModern_reduced_pruned.ind
popfilename:    OutgroupF3_French.txt
```

Now run the program qp3Pop, which comes in the AdmixTools package:

```
qp3Pop -p OutgroupF3.par > OutgroupF3_French.log
```

Look at the output log file. Which population shares the largest amount of drift with French? Which shares the least amount of drift?

We can also plot the results. We need to extract the lines from the log file that are of interest to us. We'll also select the unnecessary columns at the beginning of each line:

```
grep 'result:' OutgroupF3_French.log | awk '{print $2, $3, $4, $5, $6, $7, $8, $9}' > OutgroupF3_French.tab
```

Let's load it into R and create a plot with error bars for each of the F3 values:

```
R
library("Hmisc")
f3tab = read.table("OutgroupF3_French.tab", col.names=c("PopA", "PopB", "PopC", "F3", "SE", "Z", "SNPs"))
f3ordered = f3tab[order(-f3tab$F3),]
numrows = dim(f3ordered)[1]
errbar(1:numrows, f3ordered$F3,(f3ordered$F3+f3ordered$SE),(f3ordered$F3-f3ordered$SE), pch=20, las=2, xaxt='n',xlab="", ylab="F3")
axis(1, at=1:numrows, labels=f3ordered$PopB, las=2)
```

# Admixture F3 statistics

Now, we'll hypothesize that French are an admixed group resulting from the mixture of two populations. We'll try to find the pair of populations in our panel that can best stand in for those two populations. For this, we'll resort to using Admixture F3 statistics, in which the population of interest (French) is now placed in the third position of our population file (where Ju_hoan_North was previously when computing Outgroup F3 statistics). The first and second column will be a candidate pair of populations. The more negative the Admixture F3 statistic, the stronger the violation of "tree-ness" relating the 3 populations - an indication of admixture or population structure.

We'll use a pre-made R script ($SCRIPTS/BuildF3List.R) to iterate over all possible pairs of populations that do not include French:

```
cut -f 3 $DATA/AncientModern_reduced_pruned.ind | sort | uniq > allpops.txt
Rscript $SCRIPTS/BuildF3List.R allpops.txt French
```

We now have a new output file called "f3target_French.popfile", which always contains the French in the 3rd column and all other possible pairs of populations in the first and second columns.

Let's prepare a parameter file, this time called AdmixtureF3.par:

```
genotypename:   [YOUR DATA DIRECTORY HERE]/AncientModern_reduced_pruned.geno
snpname:	[YOUR DATA DIRECTORY HERE]/AncientModern_reduced_pruned.snp
indivname:	[YOUR DATA DIRECTORY HERE]/AncientModern_reduced_pruned.ind
popfilename:    f3target_French.popfile
```

Finally, run qp3Pop:

```
qp3Pop -p AdmixtureF3.par > AdmixtureF3_French.log
```

As before, let's extract the interesting rows and plot them in R. This time, we'll use pair lables rather than the labels of individual popualtions in the x-axis. We'll also order the populations in increasing order:

```
grep 'result:' AdmixtureF3_French.log | awk '{print $2, $3, $4, $5, $6, $7, $8, $9}' > AdmixtureF3_French.tab
R
library("Hmisc")
f3tab = read.table("AdmixtureF3_French.tab", col.names=c("PopA", "PopB", "Target", "F3", "SE", "Z", "SNPs"))
f3ordered = f3tab[order(f3tab$F3),]
pairs = paste(f3ordered$PopA,f3ordered$PopB,sep=" + ")
numrows = dim(f3ordered)[1]
errbar(1:numrows, f3ordered$F3,(f3ordered$F3+f3ordered$SE),(f3ordered$F3-f3ordered$SE), pch=20, las=2, xaxt='n',xlab="", ylab="F3")
axis(1, at=1:numrows, labels=pairs, las=2,cex.axis=0.4)
abline(h=0)
```

Are any F3 statistics negative? Are they significantly so? (|K| > 3?) Which pairs of populations do these correspond to? What could this mean?

# D-statistics

Finally, we can further confirm that there's been an admixture event in the past history of the French population by lookin for an excess of ABBA or an excess of BABA sites in a 4-population tree. Let's assume that Ju_hoan_North is a distant outgroup and that the French are a sister clade to Sardinians, possibly with some additional admixture from another group. Let's test two groups as potential sources of this admixture event: Karitiana (a Native South American group) and Yoruba (a Western African group). 

As before, we need to create a population file (Dstats_French.txt), but this time it will have to have quadruplets instead of triplets:

```
French  Sardinian       Steppe_EMBA       Ju_hoan_North
French  Sardinian       Karitiana       Ju_hoan_North
French  Sardinian       Yoruba          Ju_hoan_North
```

The parameter file (Dstats.par) should include the following information:

```
genotypename:   [YOUR DATA DIRECTORY HERE]/AncientModern_reduced_pruned.geno
snpname:        [YOUR DATA DIRECTORY HERE]/AncientModern_reduced_pruned.snp
indivname:      [YOUR DATA DIRECTORY HERE]/AncientModern_reduced_pruned.ind
popfilename:    Dstats_French.txt
```

Now, run the qpDstat program in AdmixTools:

```
qpDstat -p Dstats.par > Dstats_French.log
```

Look at the output log file. The first numerical column corresponds to the value of the D-statistic (defined here as (BABA-ABBA)/(BABA+ABBA)). The second numerical column is the Z-score corresponding to this D-statistic. The 3rd and 4th numerical columns are the BABA and ABBA counts, respectively, and the last column is the total number of SNPs available for analysis.

Are any of the two statistics significant? In what direction? (i.e. is there an excess of ABBA or BABA patterns, relative to what you would expect under a 4-population tree?). What could this mean?

To read more about this signal, you can check out the following papers, which suggest that a Northern Eurasian Steppe population that contributed ancestry to Native Americans also admixed with the ancestors of certain European populations, including French. Sardinians have almost none of this admixture, which is why F3(Sardinian,Karitiana,French) is so negative:

http://www.genetics.org/content/192/3/1065

https://www.nature.com/articles/nature12736/

https://www.nature.com/articles/nature13673

# TreeMix
Let's begin by fitting the following populations using TreeMix, an (almost-)unsupervised admixture graph program. We'll use the same populations we worked on when performing the PCA and Admixture analyses.

```
mkdir TreeMix
```

First, we need to stratify our individual allele frequencies by populations. For this, we'll use the --freq functionality in plink:

```
plink --bfile $DATA/AncientModern_reduced_pruned --freq --missing --family --out $DATA/AncientModern_reduced_pruned
gzip $DATA/AncientModern_reduced_pruned.frq.strat
```

Let's convert our plink files into treemix format.

```
python $SCRIPTS/plink2treemix.py $DATA/AncientModern_reduced_pruned.frq.strat.gz $DATA/AncientModern_reduced_pruned.treemix.frq.gz
```
Now, let's sequentially run TreeMix with 0, 1, 2 and 3 migration events. We'll set Ju_hoan_North to be on one side of the root of the tree.

```
for mig in {0,1,2,3}; do
treemix -i $DATA/AncientModern_reduced_pruned.treemix.frq.gz -o TreeMix/treemix_output_m$mig -m $mig -root Ju_hoan_North -k 50
done
```

Let's also make a file containing a list of the populations we're studying (we'll need this later to plot the residuals of our fitted graphs):

```
echo -e "Europe_LNBA\nSteppe_EMBA\nJu_hoan_North\nMbuti\nYoruba\nFrench\nSardinian\nItalian_North\nOrcadian\nPapuan\nAmi\nMayan\nKaritiana" > TreeMix/pop_order.txt
```

We can visualize the results using R scripts that can be downloaded along with the treemix program. For example, for a tree with no migration events, we can plot the corresponding graph as follows.

```
R
source("/science/groupdirs-nfs/SCIENCE-SNM-Archaeo/class_aDNA_2018/scripts/plotting_funcs.R")
plot_tree("TreeMix/treemix_output_m0")
```

Plot the other graphs and study their respective topologies. Which admixture events do you observe? Do these make sense based on your knowledge of human history? Note that some migration events may be added because of poor representation of certain populations that may have been important in human history. For example, we're not including Denisovans and Neanderthals in our graphs, which are known to have contributed ancestry to Papuans and non-Africans, respectively.

Take a look at the length of the branches in the tree. Why are some branches much longer than others? What does the length here represent?

You can also plot the residual fits from each graph. For example, for the graph containing no migrations:

```
plot_resid("TreeMix/treemix_output_m0", "TreeMix/pop_order.txt")
```

Take a look at these residuals. Which pairs of populations are worst-fitted under each graph?


