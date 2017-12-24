
# Data processing

First, let's convert our pruned plink files back to packed eigenstrat format. Use a new parameter file (plink2geno.par), with the following information:

```
genotypename: Data/HumanOriginsPublic2068_reduced_pruned.bed
snpname: Data/HumanOriginsPublic2068_reduced_pruned.bim
indivname: Data/HumanOriginsPublic2068_reduced_pruned.fam
outputformat: PACKEDANCESTRYMAP
genotypeoutname: Data/HumanOriginsPublic2068_reduced_pruned.geno
snpoutname: Data/HumanOriginsPublic2068_reduced_pruned.snp
indivoutname: Data/HumanOriginsPublic2068_reduced_pruned.ind
genotype file processed
```

Run convertf:

```
convertf -p plink2geno.par
```

Before we proceed, we need to clean up our ind file:

```

```


# Outgroup F3 statistics


# Admixture F3 statistics

# D-statistic


# qpWave / qpAdm

