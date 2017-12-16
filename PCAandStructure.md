# Record chrX and chrY SNPs
cat NearEastPublic/HumanOriginsPublic2068.snp | tr -s " " | awk 'BEGIN{OFS="\t"}{if ($2 == "23" || $2 == "24") print}' > Data/toremove.snp

mkdir Data
# Convert from packed eigenstrat to packed ped format
convertf -p geno2plink.par
# Fix *fam file
paste <(cat NearEastPublic/HumanOriginsPublic2068.ind | awk '{print $3}') <( cat Data/HumanOriginsPublic2068.fam | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$
5,$6}') > temp.fam
mv temp.fam Data/HumanOriginsPublic2068.fam
# Make list of families to extract
echo -e "Ju_hoan_North\nSardinian\nFrench\nItalian_North\tHan\nAmi\nYoruba\nMbuti\nPapuan\nOrcadian\nMayan\nKaritiana" > Data/groups_to_keep.txt
# Extract populations
plink --bfile Data/HumanOriginsPublic2068 --keep-fam Data/groups_to_keep.txt --make-bed --out Data/HumanOriginsPublic2068_reduced




# PCA - takes time
mkdir PCA
smartpca -p pca.par
# Visualize PCA
Rscript scripts/PlotPCA.R PCA/eigenvectors.txt Data/HumanOriginsPublic2068_reduced.fam PCA/PCA_World.pdf


# Admixture analysis - takes time
mkdir Admixture
cd Admixture
mkdir K3
# Run Admixture with K=3
cd K3; admixture ../../Data/HumanOriginsPublic2068_reduced.bed 3; cd ..

# Visualize results with pong
echo -e "K3_run1\t3\tK3/HumanOriginsPublic2068_reduced.3.Q" > filemap.txt
