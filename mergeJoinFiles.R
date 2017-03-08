setwd('/Users/jiemingchen/workspace/1KG_p3_LD')

## this script takes in 4 args: file1 , file2, a flag, and an output filename and does a join (merge) based on the flag
## note that colnames are hardcoded!
## flag
## 1: inner join
## 2: left outer join (file 1)
## 3: right outer join (file 2)
## 4: outer join (things not merge-able are placed NA)
## http://stackoverflow.com/questions/1299871/how-to-join-merge-data-frames-inner-outer-left-right

# args = commandArgs(trailingOnly=TRUE)

# test if number of arguments NOT 6: if not, return an error
# if (length(args)!=4) {
#   stop("Only four arguments can be supplied: file1 file2 flag outputname", call.=FALSE)
# }

# read file1 and file2
data1 = read.table("CEU-chr10-ld-hapr2-mt0.5-500mb.hap.ld.mod", header = T, sep = "\t", stringsAsFactors = FALSE, na.strings = "")
data2 = read.delim("dbsnp_138.hg19.noChr.txt", header = T, sep = "\t", stringsAsFactors = FALSE, na.strings = "")

f = 2
# merge by flag
if (f == 1)
{
  # inner join
  data12 = merge(data1, data2, by.x = args[2], by.y = args[4])
}else if (f == 2)
{
  # left outer join
  # data12 = merge(data1, data2, by.x = 'CHR.POS1', by.y = 'CHROM.POS', all.x = TRUE)
  data12_order = merge.with.order(data1, data2, by.x = 'CHR.POS1', by.y = 'CHROM.POS', all.x = TRUE, keep_order = 1)
  data123 = merge.with.order(data12_order, data2, by.x = 'CHR.POS2', by.y = 'CHROM.POS', all.x = TRUE, keep_order = 1)
}else if (f == 3)
{
  # right outer join
  data12 = merge(data1, data2, by.x = args[2], by.y = args[4], all.y = TRUE)
}else if (f == 4)
{
  # outer join
  data12 = merge(data1, data2, by.x = args[2], by.y = args[4], all = TRUE)
}

# output
write.table(data123, file="CEU-chr10-ld-hapr2-mt0.5-500mb.hap.ld.snp1.dbsnp138_", row.names=FALSE, sep = "\t", quote = FALSE)
