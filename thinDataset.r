setwd("C:/Documents and Settings/chenjm/Desktop/thinDataset")

# this file should have headers snp-id, chromosome, position
# sort the input file if it's marker file according to chr and position
data <- read.table("jm-sorted", header=T)

# this one gives number of rows in the data matrix
size=dim(data)[1]  

# k is the kth set u want to choose for the runs
# example: if i am taking every 60, I would have ~60 sets
# but I only want sets 1,20,40
# in set 1: snp 1, 61,121...
# in set 20: snp 20, 80, 140...
# in set 40: snp 40, 100, 160...
k <- c(1,3)
#k <- c(1)

for (i in k)
#for (i in 1:26)
{
	# example:
	# this samples the big set by taking every 60: 1,61,121...
	# thins into 60 sets if i=1, k=60
	# set <- seq(from=i, to=size, by=12)
	set <- seq(from=i, to=size, by=4)
	data.new <- data[set,]
	filename = paste("snpset",i, ".txt",sep="")
	write.table(data.new, file=filename, quote=F, row.names=F, sep="\t")
}
