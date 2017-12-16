
# TreeMix analysis
mkdir TreeMix
# Stratify frequencies by populations
plink --bfile Data/HumanOriginsPublic2068_reduced --freq --missing --family --out Data/HumanOriginsPublic2068_reduced
gzip Data/HumanOriginsPublic2068_reduced.frq.strat
# Convert from plink to treemix format
python scripts/plink2treemix.py Data/HumanOriginsPublic2068_reduced.frq.strat.gz Data/HumanOriginsPublic2068_reduced.treemix.frq.gz
# Run TreeMix
for mig in {0,1,2,3}; do
treemix -i Data/HumanOriginsPublic2068_reduced.treemix.frq.gz -o TreeMix/treemix_output_m$mig -m $mig -root Ju_hoan_North -k 50
done
# Visualize results


# qpGraph
