#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

table <- read.table(args[1],header=FALSE)

allpops <-  as.character(unique(table[,1]))



# Create triplets
pairs <- t(combn(allpops,2))
triplets <- c()

targets <- strsplit(args[2],",")[[1]]

for(target in targets){
addtrip <- cbind(pairs,target)
triplets <- rbind(triplets,addtrip)
}

# Remove triplets with repeats
tokeep <- apply(triplets,1,function(triplet){
return(anyDuplicated(triplet) == 0)
})
triplets <- triplets[tokeep,]

write.table(triplets,file=paste("f3target_",args[2],".popfile",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
