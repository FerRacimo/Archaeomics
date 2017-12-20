# TreeMix
Let's begin by fitting the following populations using TreeMix, an (almost-)unsupervised admixture graph program: XXX,XXX,XXX

```
mkdir TreeMix
```

First, we need to stratify our individual allele frequencies by populations:

```
plink --bfile Data/HumanOriginsPublic2068_reduced --freq --missing --family --out Data/HumanOriginsPublic2068_reduced
gzip Data/HumanOriginsPublic2068_reduced.frq.strat
```

Let's convert our plink files into treemix format.

```
python scripts/plink2treemix.py Data/HumanOriginsPublic2068_reduced.frq.strat.gz Data/HumanOriginsPublic2068_reduced.treemix.frq.gz
```
Now, let's sequentially run TreeMix with 0, 1, 2 and 3 migration events:

```
for mig in {0,1,2,3}; do
treemix -i Data/HumanOriginsPublic2068_reduced.treemix.frq.gz -o TreeMix/treemix_output_m$mig -m $mig -root Ju_hoan_North -k 50
done
```

We can visualize the results using the XXX R script.

Study the topology of the graph? Which admixture events do we observe? Do these make sense based on your knowledge of human history?

Take a look at the length of the branches in the tree. Why are some branches much longer than others? What does the length here represent?


# qpGraph

We can now do a more supervised exploration of our populations' history, using qpGraph...
