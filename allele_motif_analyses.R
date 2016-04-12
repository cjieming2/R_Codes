setwd("C:/Users/Jieming/Documents/thesis/mark_work/allele_specificity/motif_analyses")

filename = 'zintHets.asb.acc.14samples.allelicratio.motif.mode0_1.redun.txt'
data = read.table(filename, header=F, stringsAsFactors = F, comment.char="")

freq_diff = data$V8 - data$V9
data["freq_diff"] <- freq_diff
pear = cor.test(freq_diff, data$V5, method="pearson")
spear = cor.test(freq_diff, data$V5, method="spearman")


## fisher's test
q1 = sum(data$freq_diff<=0 & data$V12<=0.5)
q2 = sum(data$freq_diff>0 & data$V12<=0.5)
q3 = sum(data$freq_diff<=0 & data$V12>0.5)
q4 = sum(data$freq_diff>0 & data$V12>0.5)


freqVSaratio <- matrix(c(q1,q2,q3,q4), 2,2,
                     dimnames = list(Pathways = c("freqdiff<=0","freqdiff>0"),SNPs = c("AR<=0.5","AR>0.5")))
x = fisher.test(freqVSaratio,alternative="two.sided")
y = c("freqdiff<=0;AR<0.5",paste(x$estimate),x$p.value)