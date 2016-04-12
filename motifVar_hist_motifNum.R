setwd("C:/Users/Jieming/Documents/thesis/lynne_work/motifVar/motifSize")
source("C:/Users/Jieming/Documents/Shared/scripts-R_perl_shell_macros/R codes/jmRlib.R")

filename = 'HUMAN_TPR.prot.motifSize'
data = read.table(filename, header=F, stringsAsFactors = F, comment.char="", sep="\t")

x11()
par(cex.axis=1,cex.lab=1.4)
dtab = as.data.frame(table(data$V1))
dtab = rbind(dtab[1:7,],  
             data.frame(Var1="38",Freq="0"),
             dtab[8,])
barplot(as.integer(dtab$Freq),names.arg=dtab$Var1,
        col=c("blue"), xlab="length of TPR motif (number of residues)", ylab="number of TPR motifs")
