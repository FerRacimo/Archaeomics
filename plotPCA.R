library(optparse)
library(ggplot2)

option_list <- list(make_option(c('-a','--eigenvalues'), action='store', type='character', default=NULL, help='Eigenvalue file'),
		    make_option(c('-e','--eigenvectors'), action='store', type='character', default=NULL, help='Eigenvector file'),
                    make_option(c('-c','--comp'), action='store', type='character', default=1-2, help='Number of components'),
                    make_option(c('-f','--famfile'), action='store', type='character', default=NULL, help='fam file'),
                    make_option(c('-o','--outfile'), action='store', type='character', default=NULL, help='Output file')
                    )
opt <- parse_args(OptionParser(option_list = option_list))


# Parse components to analyze
comp <- as.numeric(strsplit(opt$comp, "-", fixed=TRUE)[[1]])

eigenvalues <- read.table(opt$eigenvalues,header=FALSE)
eigenvalues <- eigenvalues / sum(eigenvalues)
eigenvectors <- read.table(opt$eigenvectors,header=FALSE,skip=1)
eigenvectors <- eigenvectors[,-dim(eigenvectors)[2]]
pops <- read.table(opt$famfile,header=FALSE)

# Plot
PC <- as.data.frame(eigenvectors)
PC[,1] <- factor(pops[,1])
colnames(PC) <-  c("Pop",paste("PC",seq(1,dim(PC)[2]-1),sep=""))

title <- paste("PC",comp[1]," (",signif(eigenvalues[comp[1],], digits=3)*100,"%)"," / PC",comp[2]," (",signif(eigenvalues[comp[2],], digits=3)*100,"%)",sep="",collapse="")

x_axis = paste("PC",comp[1],sep="")
y_axis = paste("PC",comp[2],sep="")

pdf(opt$outfile)
ggplot() + geom_point(data=PC, aes_string(x=x_axis, y=y_axis, color="Pop")) + ggtitle(title)
dev.off()


#ggsave(opt$outfile)
#print("HERE")
#unlink("Rplots.pdf", force=TRUE)
