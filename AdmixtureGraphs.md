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

We can also do a more supervised exploration of our populations' history with qpGraph, which uses the set of all possible f4 statistics among our set of populations as summary statistics to fit a graph to the data. Let's create a folder to dump our results.

```
mkdir qpGraph
```

We need to create a parameter file (graphpar.par) first:

```
DIR:    Data
indivname:       DIR/merged.v3.fixed.ind
snpname:         DIR/merged.v3.fixed.snp
genotypename:    DIR/merged.v3.fixed.geno
outpop:         NULL
blgsize:        0.05
lsqmode:       YES
diag:          .0001
hires:         YES
initmix:      1000
precision:    .0001
zthresh:      3.0
terse:        NO
useallsnps:   NO
```

This file specifies the location of the genotype files are and some optional settings in the admixture graph fitting. You can learn more about the specific parameters in the qpGraph manual: https://github.com/DReichLab/AdmixTools/blob/master/README.QPGRAPH

We also need to create a graph topology file. Unlike in TreeMix, qpGraph requires the user to specify the topology of the graph, and finds the best fitting branch lengths and admixture rates that are possible under that particular topology. We'll focus on a subset of populations: Karitiana, Yoruba, French, Sardinian and Orcadian. Below is a (very badly fitting) example. Copy this text to a file and call it graph1.txt:

```
root    r
label   Yoruba  Yoruba
label   French  French
label   Orcadian        Orcadian
label   Sardinian       Sardinian
label   Karitiana       Karitiana
edge    Yoruba_r        r       Yoruba
edge    r_x             r       x
edge    x_b     x       b
edge    b_Sardinian     b       Sardinian
edge    x_y             x       y
edge    y_Orcadian      y       Orcadian
edge    y_a             y       a
edge    a_French        a       French
edge    k_Karitiana     k       Karitiana
admix   k               a       b       50      50
```

The first line sets the name of the root. The next few lines can serve to give alternative labels to the populations we'll analyze (in our case, we'll keep the same names as in the genotype file). Then, we specify each of the branches in the format:

edge [edge_label] [edge_parent_node] [edge_child_node]

Finally, we specify any admixed nodes we may want to introduce:

admix [child_node] [parent_1] [parent_2] [admixture_rate_1] [admixture_rate_2]

The admixture rates specified in the graph file only serve as a starting point for the fitting algorithm, and may be different in the outputted graph. Unless one has any prior information on these rates, one can just set them to be 50/50.

Finally, we run qpGraph and dump the output into a logfile:

```
qpGraph -p graphpar.par -g graph1.txt -d graph1.dot > qpGraph/logfile_graph1.txt
```

Take a look at the logfile. The outlier and worst-fitting statistics can point to possible improvements in the topology of the graph. The statistics involving 4 different populations are usually the most helpful. For example, the line:

```
Yor        Fre        Sar        Kar       0.000001    -0.022276    -0.022278     0.000690   -32.266 
```

corresponds to F4(Yoruba,French,Sardinian,Karitian). Under our model, this F4 statistic is equal to 0, but when measured on our raw data this value is negative (=-0.022276). The difference between the fitted F4 and the observed F4 is equal to -0.022278, the standard error  is 0.000690 and the Z-score for this difference is very large (=-32.266). The fact that our observed statistic is more negative than our fitted statistic means that the ABBA pattern should be stronger than the BABA pattern. We could try improving this fit by rearranging our graph to somehow make French and Sardinian closer to each other (relative to Yoruba and Karitiana), or making Yoruba and Karitiana closer to each other (relative to French and Sardinian). Given our prior knowledge of human history and our fitted TreeMix graphs, it may perhaps be more advisable to make French and Sardinian closer to each other...

To visualize the fitted graph, we first convert our outputted dot file to postscript, and then visualize it using...

```
dot -Tps graph1.dot > graph1.ps
...
```

Use your previous knowledge of population relationships obtained from TreeMix, f-statistics, Admixture and PCA to modify the graph file and look for graph topologies that provide better fits than the example (using the same set of populations). In general, a heuristic rule when fitting data is to try to find a fit in which the worst-fitting f4 statistic has a Z-score whose absolute value is less than 3. One also generally aims to find parsimonious fits: we can generate arbitrarily good fits by just filling the graph with a large number of admixture events, but such a fit may be overly complex and perhaps not very realistic. Once you have a well-fitting graph, try adding more populations.
