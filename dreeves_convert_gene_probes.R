#convert top 100 probes into genes to use in DAVID later.

setwd("/Users/gilizilberman/Dropbox/Gili/multivariate")

# Read in probe annotation
annot = read.delim(file="ALL_genes_probes.txt", header=FALSE)

# Match probes in the data set to the probe IDs
probes = colnames(dat3[,3:102])

probes2annot = match(probes, annot$V1)

# Get corresponding Gene IDs
allGeneIDs = annot$V2[probes2annot]

write.table(as.data.frame(allGeneIDs), file="ProtID.txt", sep = " ",row.names=FALSE, col.names=FALSE)

davit = read.delim(file="/Users/gilizilberman/Dropbox/Gili/multivariate/David_out.txt", header=TRUE)
