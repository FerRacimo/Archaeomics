# TreeMix
Let's begin by fitting the following populations using TreeMix, an (almost-)unsupervised admixture graph program: XXX,XXX,XXX

```
mkdir TreeMix
```

First, we need to stratify our individual allele frequencies by populations:

```
plink --bfile Data/HumanOriginsPublic2068_reduced_pruned --freq --missing --family --out Data/HumanOriginsPublic2068_reduced_pruned
gzip Data/HumanOriginsPublic2068_reduced.frq.strat
```

Let's convert our plink files into treemix format.

```
python scripts/plink2treemix.py Data/HumanOriginsPublic2068_reduced_pruned.frq.strat.gz Data/HumanOriginsPublic2068_reduced_pruned.treemix.frq.gz
```
Now, let's sequentially run TreeMix with 0, 1, 2 and 3 migration events:

```
for mig in {0,1,2,3}; do
treemix -i Data/HumanOriginsPublic2068_reduced_pruned.treemix.frq.gz -o TreeMix/treemix_output_m$mig -m $mig -root Ju_hoan_North -k 50
done
```

We can visualize the results using the XXX R script.

Study the topology of the graph. Which admixture events do you observe? Do these make sense based on your knowledge of human history?

Take a look at the length of the branches in the tree. Why are some branches much longer than others? What does the length here represent?


# qpGraph

We can now do a more supervised exploration of our populations' history, using qpGraph...

We need to create a parameter file first:


```
```

We also need to create a graph topology file. Here's a (badly fitting) example:

```
```

Modify the graph file until you find a topology that gives you a good fit. A heuristic rule is to find a fit in which the worst-fitting f4 statistic has a Z-score whose absolute value is less than 3. One also generally aims to find parsimonious fits: we can generate arbitrarily good fits by just filling the graph with a large number of admixture events, but such a fit may be overly complex and perhaps not very realistic.
