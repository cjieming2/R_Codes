setwd("C:/Documents and Settings/chenjm/Desktop/")

# a.txt file should have only the column for the values
# this script cannot remove zeroes, remove them in UNIX b4 doing this stats summ
# type 8 quantile calc is unbiased on any distribution of the data
# requires a header

filename = "hetero.txt"
samples <- read.table(filename, header=TRUE)

norm.interval = function(data, variance = var(data), conf.level = 0.95) {
z = qnorm((1 - conf.level)/2, lower.tail = FALSE)
xbar = mean(data)
sdx = sqrt(variance/length(t(data)))
c(xbar - z * sdx, xbar + z * sdx)
}

summary.stats = function(data) {
cat("<SUMMARY STATISTICS> for", filename, "\n")
cat("Minimum value: ", min(data), "\n")
cat("Maximum value: ", max(data), "\n")
cat("\n---QUANTILES---\n")
cat("10%\t25%\t50%\t75%\n")
cat(quantile(data, probs=c(10,25,50,75)/100, type=8, names = TRUE), "\n")
cat("Median: ", median(data), "\n")
cat("\n---MEAN & 95% CONFIDENCE INTERVAL OF MEAN---\n")
cat("Mean: ", mean(data), "\n")
cat("Variance: ", var(data), "\n")
cat("Standard deviation: ", sd(data), "\n")
cat("95% Confidence Interval: (", norm.interval(data), ")\n")
}

#summary.stats(samples$sample.call.rate)
summary.stats(samples$hetero)