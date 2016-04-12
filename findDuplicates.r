setwd("C:/Documents and Settings/chenjm/Desktop")

data <- read.table("900Ysnps-affy6-sgvp.mk", header=T)

## this finds the duplicates within the column 3 of the data, 
## duplicates are compared row-wise, i.e whether line 1 is a duplicate of line 2
## in field (col) 3


duplicated(data[c(10)], fromLast=TRUE)