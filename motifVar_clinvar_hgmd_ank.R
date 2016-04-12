# setwd("C:/Users/Jieming/Documents/thesis/lynne_work/motifVar/tpr/clinvar_hgmd")
setwd("C:/Users/Jieming/Documents/thesis/lynne_work/motifVar/ank/clinvar_hgmd")
source("C:/Users/Jieming/Documents/Shared/scripts-R_perl_shell_macros/R codes/jmRlib.R")

filename = 'clinvar_hgmd_analyses_anks v1.txt'
data = read.table(filename, header=T, stringsAsFactors = F, comment.char="", sep="\t")

x11()
par(cex.axis=1,cex.lab=1.4)
barplot(rbind(data$clinvar,data$hgmd.2015,data$combined.2015),names.arg=data$sorted,beside=T,
        col=c("blue","lightblue","navy"), legend=c("ClinVar","HGMD","combined"), args.legend=list(x=120),
        xlab="positions (ordered)", ylab="number of SNVs")
