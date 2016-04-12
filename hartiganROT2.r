calcHartiganROT <- function(ssK, ssKplus1, n, k){
  ((ssK/ssKplus1)-1)*(n-k-1)
}

rowCount <- nrow([DATA])
hartiganROT <- c()
comparisons <- c()
ss <- c()
for(k in 2:rowCount-1) {
  clusters <- kmeans([DATA], centers=k) 
  ssK <- clusters$tot.withinss
  ss <- c(ss, ssK)
  if(k > 2) {
    rot <- calcHartiganROT(ssKMinus1, ssK, numGenes, k)
    hartiganROT <- c(hartiganROT, rot)
    comparisons <- c(comparisons, paste(k-1, 'v', k))
  }
  ssKMinus1 <- ssK
}
hartiganRatios <- data.frame(hartiganROT=hartiganROT, row.names=comparisons)