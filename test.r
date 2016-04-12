#nchar - counts num of char

groups <- seq(from=0, to=1, by=0.05)
chin.cat <- cut(maf, breaks=groups, include.lowest=T, right=T)
summary(chin.cat)


#barplot

rounding down
x[i] = trunc((trunc(x[i]/10^((nchar(x[i]))-1)) * (10^((nchar(x[i]))-1)))/1000)

Get percentile
Fn(5000) = 0.4372