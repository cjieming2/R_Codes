#################################################################
## this script is adapted from Matthew Kan's script in plotting 
## Patient-overlap in AD-studies

setwd("/Users/jiemingchen/Documents/transplantation/a_donor/immport")

library(RImmPort) 
library(DBI) 
library(sqldf) 
library(plyr)
library(RMySQL)
library(dplyr)

## input data
mydata = read.table("studies_overlapping_samples.txt", header = T, sep = "\t", stringsAsFactors = FALSE)

## preparing the dataset
studies = as.data.frame(unique(mydata$STUDYID), stringsAsFactors = FALSE ); names(studies) = "id"
studies$size = apply(studies, 1, function(x) nrow(mydata[mydata$STUDYID == x,]))
studies_subj = lapply(studies$id, function(x) mydata[mydata$STUDYID == x,]$USUBJID)
names(studies_subj) = studies$id

## form overlapping matrix of overlapping number of individuals
a = lapply(studies_subj, function(x) lapply(studies_subj, function(y) length(intersect(x,y))))
b = do.call(cbind, a) ## note that this is a 2D matrix

a_p = lapply(studies_subj, function(x) lapply(studies_subj, function(y) length(intersect(x,y))/length(x)*100 ))
b_p = do.call(cbind, a_p) ## note that this is a 2D matrix

##############################
## plot
## Square Pie function
squarePie <- function(pct, col="black", col.grid="#e0e0e0", col.border="black", main="") {
  
  if (pct > 100) {
    pct <- 100
    warning("Percentage value, pct, should be an integer between 0 and 100")
  } else if (pct < 0) {
    pct <- 0
    warning("Percentage value, pct, should be an integer between 0 and 100.")
  }
  
  # Round to nearest integer
  pct <- round(pct)
  
  # x- and y-coordinates of rows and columns
  x_row <- 1:10
  y_col <- 1:10
  
  # put together full coordinate vectors
  x <- rep(x_row, 10)
  y <- rep(y_col, each=10)
  
  # set colors
  fill_col <- c(rep(col, pct), rep("#ffffff", 100-pct))
  
  # plot
  plot(0, 0, type="n", xlab="", ylab="", main=main, xlim=c(0, 11), ylim=c(0, 10.5), asp=1, bty="n", axes=FALSE)
  symbols (x, y, asp=1, squares=rep(1,100), inches=FALSE, add=TRUE, bg=fill_col, fg=col.grid, lwd=0.5)
  rect(.5, .5, 10.5, 10.5, lwd=2, border=col.border)
}

## Generating a matrix of square Pis
x11(type="cairo")
par(mfrow = c(nrow(studies), nrow(studies)), mar = c(0, 0, 0, 0))

mapply(squarePie, b_p, col.grid=NA, col="#007CBE", col.border="#B4B9BF")


## original code cant weave the coloring of diag vs nondiag into mapply
## for later improvements?
# for(i in 1:10) {
#   for(j in 1:10) {
#     if (i!=j) {
#       border_color = "#B4B9BF"
#       squarePie(percentage_overlap[i, j], col.grid=NA, col="#007CBE", col.border=border_color)
#     } else {
#       border_color = "#000000"
#       squarePie(percentage_overlap[i, j], col.grid=NA, col="#e0e0e0", col.border=border_color)
#     }
#     if (overlaps[i, j] > 0) {
#       text(5.5, 5.5, overlaps[i, j], font=2, cex=1.5, family = "Roboto")
#     }
#   }
# }


