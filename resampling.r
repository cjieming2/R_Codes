setwd("C:/Documents and Settings/chenjm/Desktop")

####this does resampling via jackknife or bootstrapping
## coefficent of variation
CV <- function(data) 
{
	sqrt(var(data))/mean(data)
}

## read in data
data <- read.table("sampling-test.txt", header=TRUE, sep="\t")
firstcol <- c(data)

