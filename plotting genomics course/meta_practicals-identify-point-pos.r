### Reading in the output from individual case-control analysis ### 
ceu.in <- read.table("D:/GGCP_April2010/Practicals/Day4_Session3/practical_meta_CEU_BRCA2.txt", sep="\t", header = T)
asi.in <- read.table("D:/GGCP_April2010/Practicals/Day4_Session3/practical_meta_JPT+CHB_BRCA2.txt", sep="\t", header = T)
yri.in <- read.table("D:/GGCP_April2010/Practicals/Day4_Session3/practical_meta_YRI_BRCA2.txt", sep="\t", header = T)

### Explore one of the datasets, say for the Europeans ###
names(ceu.in)
ceu.in[1:10,]

### Plot the signal plot of the genomic region in Europeans ###
plot(ceu.in[,"pos"]/10^6, -log(ceu.in[,"frequentist_add"], base = 10), bty = "n", xlab = "Physical position (Mb)", ylab = "-log10 P-value", type = "n", main = "Europeans", las = 1)
points(ceu.in[,"pos"]/10^6, -log(ceu.in[,"frequentist_add"], base = 10), cex = 0.8, col = 2, pch = 16)
identify(ceu.in[,"pos"]/10^6, -log(ceu.in[,"frequentist_add"], base = 10), labels = as.character(ceu.in[,"rsid"]))

# East Asians #
plot(asi.in[,"pos"]/10^6, -log(asi.in[,"frequentist_add"], base = 10), bty = "n", xlab = "Physical position (Mb)", ylab = "-log10 P-value", type = "n", main = "East Asians", las = 1)
points(asi.in[,"pos"]/10^6, -log(asi.in[,"frequentist_add"], base = 10), cex = 0.8, col = 3, pch = 16)
identify(asi.in[,"pos"]/10^6, -log(asi.in[,"frequentist_add"], base = 10), labels = as.character(asi.in[,"rsid"]))

# Yoruba Africans #
plot(yri.in[,"pos"]/10^6, -log(yri.in[,"frequentist_add"], base = 10), bty = "n", xlab = "Physical position (Mb)", ylab = "-log10 P-value", type = "n", main = "Yoruba Africans", las = 1)
points(yri.in[,"pos"]/10^6, -log(yri.in[,"frequentist_add"], base = 10), cex = 0.8, col = 4, pch = 16)
identify(yri.in[,"pos"]/10^6, -log(yri.in[,"frequentist_add"], base = 10), labels = as.character(yri.in[,"rsid"]))

### Identify the common SNPs in Europeans and East Asians analyses ###
common.snps <- intersect(ceu.in[,"rsid"], asi.in[,"rsid"])
dim(ceu.in)[1] # number of SNPs from the European scan 
dim(asi.in)[1] # number of SNPs from the East Asian scan
length(common.snps) # number of common SNPs
ceu.meta <- ceu.in[which(is.element(ceu.in[,"rsid"], common.snps)),]
asi.meta <- asi.in[which(is.element(asi.in[,"rsid"], common.snps)),]

# Programming function for fixed-effects inverse variance meta-analysis # 
integrand <- function(x) {x^(-0.5) * exp(-x/2) / sqrt(2 * pi)}
OR.merge <- function(input){
   K <- dim(input)[1]
   psi <- log(input[,1])
   sigma_sq <- ((psi - log(input[,2]) )/1.96)^2
   weight_sq <- 1/sigma_sq
   Psi <- sum(psi * weight_sq)/sum(weight_sq)
   var.Psi <- 1/sum(weight_sq)
   test.stat <- Psi^2/var.Psi
   Combined.OR.CI <- exp(c(Psi - 1.96 * sqrt(var.Psi), Psi + 1.96 * sqrt(var.Psi)))
   p.value <- 1 - pchisq(test.stat, 1)
   output <- list(Chisq = test.stat, p.value = p.value, OR = exp(Psi), OR_CI = Combined.OR.CI)
   output
}
heterogeneity.OR <- function(input){
   integrand <- function(x) {x^(-0.5) * exp(-x/2) / sqrt(2 * pi)}
   K <- dim(input)[1]
   psi <- log(input[,1])
   sigma_sq <- ((psi - log(input[,2]) )/1.96)^2
   weight_sq <- 1/sigma_sq
   Psi <- sum(psi * weight_sq)/sum(weight_sq)
   Q <- sum(weight_sq * (psi - Psi)^2)
   output <- list(Chisq = Q, p.value = 1 - pchisq(Q, K - 1))
   output
}
########################################################


### Perform the meta-analysis between 2 populations ###
pop1.meta <- ceu.meta
pop2.meta <- asi.meta

n.snp <- dim(pop1.meta)[1]
output <- as.data.frame(matrix(, nr = n.snp, nc = 20))
colnames(output) <- c("rsID", "position", "alleleA", "alleleB", "pop1_controls_maf", "pop1_cases_maf", "pop1_OR_alleleB", "pop1_OR_alleleB_lower", "pop1_OR_alleleB_upper", "pop1_p", "pop2_controls_maf", "pop2_cases_maf", "pop2_OR_alleleB", "pop2_OR_alleleB_lower", "pop2_OR_alleleB_upper", "pop2_p", "meta_OR_alleleB", "meta_OR_alleleB_lower", "meta_OR_alleleB_upper", "meta_p")
output[,1:4] <- pop1.meta[,1:4]
output[,5:10] <- pop1.meta[,13:18]
output[,11:16] <- pop2.meta[,13:18]
for (i in 1:n.snp){
   temp.in <- matrix(c(unlist(pop1.meta[i,c("all_OR", "all_OR_lower", "all_OR_upper")]), unlist(pop2.meta[i,c("all_OR", "all_OR_lower", "all_OR_upper")])), nr = 2, nc = 3, byrow = T)
   temp.out <- OR.merge(temp.in)
   output[i,17] <- temp.out$OR
   output[i,18:19] <- temp.out$OR_CI
   output[i,20] <- temp.out$p.value
   if (output[i,20] == 0 & is.na(output[i,20]) == F) output[i,20] <- integrand(temp.out$Chisq)
   if (is.na(output[i,20]) == T) output[i,20] <- -1   
}
rm.flag <- which(output[,20] == -1)
if (length(rm.flag) > 0) output <- output[-rm.flag,]


### Plot the combined signals in the same plot ###
plot(output[,"position"]/10^6, -log(output[,"meta_p"], base = 10), ylim = c(0, max(-log(output[,c("pop1_p","pop2_p","meta_p")], base = 10), na.rm = T)), bty = "n", xlab = "Physical position (Mb)", ylab = "-log10 P-value", type = "n", main = "Meta-analysis", las = 1, cex.main = 1.2, cex.lab = 1.2)
for (i in 1:n.snp){
   points(rep(output[i,"position"]/10^6, 3), -log(unlist(output[i,c("pop1_p", "pop2_p", "meta_p")]), base = 10), cex = 0.8, col = c(2,3,1), pch = 16, type = "b", lty = 3)
}

# Zoom into a specific region #
start.region = 31.7 
end.region = 32

flag <- which(output[,"position"]/10^6 > start.region & output[,"position"]/10^6 < end.region)
plot(output[flag,"position"]/10^6, -log(output[flag,"meta_p"], base = 10), ylim = c(0, max(-log(output[,c("pop1_p","pop2_p","meta_p")], base = 10), na.rm = T)), bty = "n", xlab = "Physical position (Mb)", ylab = "-log10 P-value", type = "n", main = "Meta-analysis", las = 1, cex.lab = 1.2, cex.main = 1.2)
for (i in 1:length(flag)){
   points(rep(output[flag[i],"position"]/10^6, 3), -log(unlist(output[flag[i],c("pop1_p", "pop2_p", "meta_p")]), base = 10), cex = 0.8, col = c(2,3,1), pch = 16, type = "b", lty = 3)
}
legend(locator(1), c("Europeans only", "East Asians only", "Meta-analysis"), pch = 16, col = c(2,3,1), bty = "n", cex = 1.2)
identify(output[flag,"position"]/10^6, -log(output[flag,"meta_p"], base = 10), labels = as.character(output[flag, "rsID"]))


### Perform all possible meta-analyses between 3 populations ###
common.snps <- intersect(intersect(ceu.in[,"rsid"], asi.in[,"rsid"]), yri.in[,"rsid"])
dim(ceu.in)[1] # number of SNPs from the European scan 
dim(asi.in)[1] # number of SNPs from the East Asian scan
dim(yri.in)[1] # number of SNPs from the Yoruba African scan
length(common.snps) # number of common SNPs
ceu.meta <- ceu.in[which(is.element(ceu.in[,"rsid"], common.snps)),]
asi.meta <- asi.in[which(is.element(asi.in[,"rsid"], common.snps)),]
yri.meta <- yri.in[which(is.element(yri.in[,"rsid"], common.snps)),]

pop1.meta <- ceu.meta
pop2.meta <- asi.meta
pop3.meta <- yri.meta

n.snp <- dim(pop1.meta)[1]
output <- as.data.frame(matrix(, nr = n.snp, nc = 38))
colnames(output) <- c("rsID", "position", "alleleA", "alleleB", "pop1_controls_maf", "pop1_cases_maf", "pop1_OR_alleleB", "pop1_OR_alleleB_lower", "pop1_OR_alleleB_upper", "pop1_p", "pop2_controls_maf", "pop2_cases_maf", "pop2_OR_alleleB", "pop2_OR_alleleB_lower", "pop2_OR_alleleB_upper", "pop2_p", "pop3_controls_maf", "pop3_cases_maf", "pop3_OR_alleleB", "pop3_OR_alleleB_lower", "pop3_OR_alleleB_upper", "pop3_p", "meta12_OR_alleleB", "meta12_OR_alleleB_lower", "meta12_OR_alleleB_upper", "meta12_p", "meta13_OR_alleleB", "meta13_OR_alleleB_lower", "meta13_OR_alleleB_upper", "meta13_p", "meta23_OR_alleleB", "meta23_OR_alleleB_lower", "meta23_OR_alleleB_upper", "meta23_p", "meta123_OR_alleleB", "meta123_OR_alleleB_lower", "meta123_OR_alleleB_upper", "meta123_p")
output[,1:4] <- pop1.meta[,1:4]
output[,5:10] <- pop1.meta[,13:18]
output[,11:16] <- pop2.meta[,13:18]
output[,17:22] <- pop3.meta[,13:18]
for (i in 1:n.snp){
   temp.in <- matrix(c(unlist(pop1.meta[i,c("all_OR", "all_OR_lower", "all_OR_upper")]), unlist(pop2.meta[i,c("all_OR", "all_OR_lower", "all_OR_upper")])), nr = 2, nc = 3, byrow = T)
   temp.out <- OR.merge(temp.in)
   output[i,23] <- temp.out$OR
   output[i,24:25] <- temp.out$OR_CI
   output[i,26] <- temp.out$p.value
   if (output[i,26] == 0 & is.na(output[i,26]) == F) output[i,26] <- integrand(temp.out$Chisq)
   if (is.na(output[i,26]) == T) output[i,26] <- -1   
   temp.in <- matrix(c(unlist(pop1.meta[i,c("all_OR", "all_OR_lower", "all_OR_upper")]), unlist(pop3.meta[i,c("all_OR", "all_OR_lower", "all_OR_upper")])), nr = 2, nc = 3, byrow = T)
   temp.out <- OR.merge(temp.in)
   output[i,27] <- temp.out$OR
   output[i,28:29] <- temp.out$OR_CI
   output[i,30] <- temp.out$p.value
   if (output[i,30] == 0 & is.na(output[i,30]) == F) output[i,30] <- integrand(temp.out$Chisq)
   if (is.na(output[i,30]) == T) output[i,30] <- -1   
   temp.in <- matrix(c(unlist(pop2.meta[i,c("all_OR", "all_OR_lower", "all_OR_upper")]), unlist(pop3.meta[i,c("all_OR", "all_OR_lower", "all_OR_upper")])), nr = 2, nc = 3, byrow = T)
   temp.out <- OR.merge(temp.in)
   output[i,31] <- temp.out$OR
   output[i,32:33] <- temp.out$OR_CI
   output[i,34] <- temp.out$p.value
   if (output[i,34] == 0 & is.na(output[i,34]) == F) output[i,34] <- integrand(temp.out$Chisq)
   if (is.na(output[i,34]) == T) output[i,34] <- -1   
   temp.in <- matrix(c(unlist(pop1.meta[i,c("all_OR", "all_OR_lower", "all_OR_upper")]), unlist(pop2.meta[i,c("all_OR", "all_OR_lower", "all_OR_upper")]), unlist(pop3.meta[i,c("all_OR", "all_OR_lower", "all_OR_upper")])), nr = 3, nc = 3, byrow = T)
   temp.out <- OR.merge(temp.in)
   output[i,35] <- temp.out$OR
   output[i,36:37] <- temp.out$OR_CI
   output[i,38] <- temp.out$p.value
   if (output[i,38] == 0 & is.na(output[i,38]) == F) output[i,38] <- integrand(temp.out$Chisq)
   if (is.na(output[i,38]) == T) output[i,38] <- -1   
}
rm.flag <- which(output[,38] == -1)
if (length(rm.flag) > 0) output <- output[-rm.flag,]

# Plot the signal plot #
plot(output[,"position"]/10^6, -log(output[,"meta123_p"], base = 10), ylim = c(0, max(-log(output[,c("pop1_p","pop2_p","pop3_p","meta123_p")], base = 10), na.rm = T)), bty = "n", xlab = "Physical position (Mb)", ylab = "-log10 P-value", type = "n", main = "Meta-analysis", las = 1, cex.lab = 1.2, cex.main = 1.2)
for (i in 1:n.snp){
   points(rep(output[i,"position"]/10^6, 4), -log(unlist(output[i,c("pop1_p", "pop2_p", "pop3_p", "meta123_p")]), base = 10), cex = 0.8, col = c(2,3,4,1), pch = 16, type = "b", lty = 3)
}

# Zoom into a specific region #
start.region = 31.7 
end.region = 32

flag <- which(output[,"position"]/10^6 > start.region & output[,"position"]/10^6 < end.region)
plot(output[flag,"position"]/10^6, -log(output[flag,"meta123_p"], base = 10), ylim = c(0, max(-log(output[,c("pop1_p","pop2_p","pop3_p","meta123_p")], base = 10), na.rm = T)), bty = "n", xlab = "Physical position (Mb)", ylab = "-log10 P-value", type = "n", main = "Meta-analysis", las = 1, cex.lab = 1.2, cex.main = 1.2)
for (i in 1:length(flag)){
   points(rep(output[flag[i],"position"]/10^6, 4), -log(unlist(output[flag[i],c("pop1_p", "pop2_p", "pop3_p", "meta123_p")]), base = 10), cex = 0.8, col = c(2,3,4,1), pch = 16, type = "b", lty = 3)
}
legend(locator(1), c("Europeans only", "East Asians only", "Yoruba Africans only", "Meta-analysis"), pch = 16, col = c(2,3,4,1), bty = "n", cex = 1.2)
identify(output[flag,"position"]/10^6, -log(output[flag,"meta123_p"], base = 10), labels = as.character(output[flag, "rsID"]))

# Produce a forest plot of the effect sizes for a specific SNP #
index.snp = "rs9634672"

index.flag <- which(output[,"rsID"] == index.snp)
plot(0, 1, type = "n", xlim = c(min(0.5, min(output[index.flag,c("pop1_OR_alleleB_lower", "pop2_OR_alleleB_lower", "pop3_OR_alleleB_lower", "meta12_OR_alleleB_lower", "meta13_OR_alleleB_lower", "meta23_OR_alleleB_lower", "meta123_OR_alleleB_lower")], na.rm = T)), ceiling(max(output[index.flag,c("pop1_OR_alleleB_upper", "pop2_OR_alleleB_upper", "pop3_OR_alleleB_upper", "meta12_OR_alleleB_upper", "meta13_OR_alleleB_upper", "meta23_OR_alleleB_upper", "meta123_OR_alleleB_upper")], na.rm = T))), ylim = c(0,7), yaxt = "n", bty = "n", xlab = "Odds ratio", ylab = "Population", cex.lab = 1.2, main = paste("Forest plot for ", index.snp, sep=""), cex.main = 1.2)
axis(side = 2, at = 7:1, labels = c("1", "2", "3", "1-2", "1-3", "2-3", "1-2-3"), las = 1, cex.lab = 1.2)
abline(v = 1, lty = 3)
points(output[index.flag, c("pop1_OR_alleleB_lower", "pop1_OR_alleleB_upper")], c(7,7), type = "l", lwd = 2)
points(output[index.flag, "pop1_OR_alleleB"], 7, pch = 16, cex = 1.5)
points(output[index.flag, c("pop2_OR_alleleB_lower", "pop2_OR_alleleB_upper")], c(6,6), type = "l", lwd = 2)
points(output[index.flag, "pop2_OR_alleleB"], 6, pch = 16, cex = 1.5)
points(output[index.flag, c("pop3_OR_alleleB_lower", "pop3_OR_alleleB_upper")], c(5,5), type = "l", lwd = 2)
points(output[index.flag, "pop3_OR_alleleB"], 5, pch = 16, cex = 1.5)
points(output[index.flag, c("meta12_OR_alleleB_lower", "meta12_OR_alleleB_upper")], c(4,4), type = "l", lwd = 2)
points(output[index.flag, "meta12_OR_alleleB"], 4, pch = 18, cex = 2)
points(output[index.flag, c("meta13_OR_alleleB_lower", "meta13_OR_alleleB_upper")], c(3,3), type = "l", lwd = 2)
points(output[index.flag, "meta13_OR_alleleB"], 3, pch = 18, cex = 2)
points(output[index.flag, c("meta23_OR_alleleB_lower", "meta23_OR_alleleB_upper")], c(2,2), type = "l", lwd = 2)
points(output[index.flag, "meta23_OR_alleleB"], 2, pch = 18, cex = 2)
points(output[index.flag, c("meta123_OR_alleleB_lower", "meta123_OR_alleleB_upper")], c(1,1), type = "l", lwd = 2)
points(output[index.flag, "meta123_OR_alleleB"], 1, pch = 18, cex = 3)

##########
ceu.in <- read.table("D:/GGCP_April2010/Practicals/Day4_Session3/practical_meta_CEU_CDKAL1.txt", sep="\t", header = T)
asi.in <- read.table("D:/GGCP_April2010/Practicals/Day4_Session3/practical_meta_JPT+CHB_CDKAL1.txt", sep="\t", header = T)
yri.in <- read.table("D:/GGCP_April2010/Practicals/Day4_Session3/practical_meta_YRI_CDKAL1.txt", sep="\t", header = T)

ceu.in <- read.table("D:/GGCP_April2010/Practicals/Day4_Session3/practical_meta_CEU_chr12.txt", sep="\t", header = T)
asi.in <- read.table("D:/GGCP_April2010/Practicals/Day4_Session3/practical_meta_JPT+CHB_chr12.txt", sep="\t", header = T)
yri.in <- read.table("D:/GGCP_April2010/Practicals/Day4_Session3/practical_meta_YRI_chr12.txt", sep="\t", header = T)


### Heterogeneity tests ### 
index.snp = "rs9634672"

index.flag <- which(output[,"rsID"] == index.snp)
pop.effect.sizes <- matrix(c(unlist(pop1.meta[index.flag,c("all_OR", "all_OR_lower", "all_OR_upper")]), unlist(pop2.meta[index.flag,c("all_OR", "all_OR_lower", "all_OR_upper")]), unlist(pop3.meta[index.flag,c("all_OR", "all_OR_lower", "all_OR_upper")])), nr = 3, nc = 3, byrow = T)

heterogeneity.OR(pop.effect.sizes)

