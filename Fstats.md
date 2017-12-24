
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

We'll begin by computing outgroup F3 statistics to determine which are the populations that share the most drift with French. We'll use Ju_hoan_North as the outgroup.

# Admixture F3 statistics

Now, let's hypothesize that French are an admixed group. We'll try to find...

# D-statistics


