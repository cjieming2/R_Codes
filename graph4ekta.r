require(igraph)
# setwd("C:/Users/JM/Desktop")
setwd("C:/Users/JM/thesis/mark_work/ekta_netsnps/correlations_hgnc_0deg")

## input
file = "pairs.txt"
data <- as.matrix(read.table(file, header=TRUE, sep="\t"))

mygraph <- graph.edgelist(data)

mydeg <- data.frame(cbind(degree(mygraph,mode="total"),degree(mygraph,mode="in"),degree(mygraph,mode="out")),
                    row.names=V(mygraph)$name)
names(mydeg) <- c("total","in","out")
write.table(mydeg,file=paste(file,"mygraph.txt",sep=""),col.names=NA)
