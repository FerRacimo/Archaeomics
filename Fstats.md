
# Data processing

First, let's convert our pruned plink files to eigenstrat format. Use a new parameter file (plink2geno.par), with the following information:

```
genotypename: Data/HumanOriginsPublic2068_reduced_pruned.bed
snpname: Data/HumanOriginsPublic2068_reduced_pruned.bim
indivname: Data/HumanOriginsPublic2068_reduced_pruned.fam
outputformat: EIGENSTRAT
genotypeoutname: Data/HumanOriginsPublic2068_reduced_pruned.geno
snpoutname: Data/HumanOriginsPublic2068_reduced_pruned.snp
indivoutname: Data/HumanOriginsPublic2068_reduced_pruned.ind
```

Run convertf:

```
convertf -p plink2geno.par
```

Before we proceed, we need to clean up our ind file:

```
paste <(cut -d ":" -f 2 Data/HumanOriginsPublic2068_reduced_pruned.ind | cut -d " "  -f 1) <(cut -d ":" -f 2 Data/HumanOriginsPublic2068_reduced_pruned.ind | cut -d " " -f 2) <(cut -d ":" -f 1  Data/HumanOriginsPublic2068_reduced_pruned.ind | sed 's/ //g') > Data/temp.txt
mv Data/temp.txt Data/HumanOriginsPublic2068_reduced_pruned.ind
```


# Outgroup F3 statistics

First, create a folder to place our output files:

```
mkdir OutgroupF3
```

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
genotypename:   Data/HumanOriginsPublic2068_reduced_pruned.geno
snpname:	Data/HumanOriginsPublic2068_reduced_pruned.snp
indivname:	Data/HumanOriginsPublic2068_reduced_pruned.ind
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
library("Hmisc")
f3tab = read.table("OutgroupF3_French.tab", col.names=c("PopA", "PopB", "PopC", "F3", "SE", "Z", "SNPs"))
f3ordered = f3tab[order(-f3tab$F3),]
numrows = dim(f3ordered)[2]
errbar(1:numrows, f3ordered$F3,(f3ordered$F3+f3ordered$StdErr),(f3ordered$F3-f3ordered$StdErr), pch=20, las=2, cex.axis=0.4, xaxt='n',xlab="population", ylab="F3")
axis(1, at=1:numrows, labels=f3ordered$PopB, las=2, cex.axis=0.6)
```

# Admixture F3 statistics

Now, let's hypothesize that French are an admixed group. We'll try to find...

# D-statistics


