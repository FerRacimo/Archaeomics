# TreeMix
Let's begin by fitting the following populations using TreeMix, an (almost-)unsupervised admixture graph program. We'll use the same populations we worked on when performing the PCA and Admixture analyses.

```
mkdir TreeMix
```

First, we need to stratify our individual allele frequencies by populations. For this, we'll use the --freq functionality in plink:

```
plink --bfile Data/HumanOriginsPublic2068_reduced_pruned --freq --missing --family --out Data/HumanOriginsPublic2068_reduced_pruned
gzip Data/HumanOriginsPublic2068_reduced_pruned.frq.strat
```

Let's convert our plink files into treemix format.

```
python scripts/plink2treemix.py Data/HumanOriginsPublic2068_reduced_pruned.frq.strat.gz Data/HumanOriginsPublic2068_reduced_pruned.treemix.frq.gz
```
Now, let's sequentially run TreeMix with 0, 1, 2 and 3 migration events. We'll set Ju_hoan_North to be on one side of the root of the tree.

```
for mig in {0,1,2,3}; do
treemix -i Data/HumanOriginsPublic2068_reduced_pruned.treemix.frq.gz -o TreeMix/treemix_output_m$mig -m $mig -root Ju_hoan_North -k 50
done
```

Let's also make a file containing a list of the populations we're studying (we'll need this later to plot the residuals of our fitted graphs):

```
echo -e "Ju_hoan_North\nMbuti\nYoruba\nFrench\nSardinian\nItalian_North\nOrcadian\nPapuan\nAmi\nMayan\nKaritiana" > TreeMix/pop_order.txt
```

We can visualize the results using R scripts that can be downloaded along with the treemix program. For example, for a tree with no migration events, we can plot the corresponding graph as follows.

```
R
source("scripts/plotting_funcs.R")
plot_tree("TreeMix/treemix_output_m0")
```

Plot the other graphs and study their respective topologies. Which admixture events do you observe? Do these make sense based on your knowledge of human history? Note that some migration events may be added because of poor representation of certain populations that may have been important in human history. For example, we're not including Denisovans and Neanderthals in our graphs, which are known to have contributed ancestry to Papuans and non-Africans, respectively.

Take a look at the length of the branches in the tree. Why are some branches much longer than others? What does the length here represent?

You can also plot the residual fits from each graph. For example, for the graph containing no migrations:

```
R
source("scripts/plotting_funcs.R")
plot_resid("TreeMix/treemix_output_m0", "TreeMix/pop_order.txt")
```

Take a look at these residuals. Which pairs of populations are worst-fitted under each graph?


# qpGraph

We can also do a more supervised exploration of our populations' history, using qpGraph...

```
mkdir qpGraph
```

We need to create a parameter file (graphpar.par) first:

```
```

We also need to create a graph topology file. Below is a (badly fitting) example. Copy this text to a file and call it graph1.txt:

```
```

Finally, we run qpGraph and dump the output into a logfile:

```
qpGraph -p graphpar.par -g graph1.txt > qpGraph/logfile_graph1.txt
```

Take a look at the logfile...

Modify the graph file and find topologies that provide better fits than the example. In general, a heuristic rule when fitting data is to try to find a fit in which the worst-fitting f4 statistic has a Z-score whose absolute value is less than 3. One also generally aims to find parsimonious fits: we can generate arbitrarily good fits by just filling the graph with a large number of admixture events, but such a fit may be overly complex and perhaps not very realistic.
