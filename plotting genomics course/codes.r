## ========================================================================== ##
## Genetics and Genomics of Complex Phenotypes April 2010
## Day3_Session5 (Analysing a GWAS dataset)
## ========================================================================== ##

## Section 1 ##
## --------- ##
# 1a. Reading in the data #
dat <- read.table("D:/GGCP_April2010/Practicals/Day3_Session4/case_control_10snps.txt", header=T, sep="\t", stringsAsFactors=F)
head(dat) # Explore the data.
dim(dat)

# Cross-tabulations #
snp <- "rs1080013"   # Change this to look at other SNPs.

snp.col <- which(colnames(dat)%in%snp)
ftable(dat[,"affection"], dat[,snp.col], dnn=c("affection",snp))

# 1b. Chisq-test #
chisq.test(dat[,"affection"], dat[,snp.col])

# 1c. Logistic regression -- additive model #
model1 <- glm(dat[,"affection"] ~ dat[,snp.col], family="binomial")  # Fit the model
summary(model1) # look at the model summary statistics

# 1d. Logistic regression -- general model #
model2 <- glm(dat[,"affection"] ~ as.factor(dat[,snp.col]), family="binomial")  # Fit the model
summary(model2) # look at the model summary statistics
1-pchisq(abs(model2$deviance - model2$null.deviance), df=abs(model2$df.null - model2$df.residual))

# 1e. Logistic regression -- additive model, adjusting for gender #
model3 <- glm(dat[,"affection"] ~ dat[,snp.col] + dat[,"gender"]-1, family="binomial")  # Fit the model
summary(model3)

# 1f. Looping over all the SNPs, fitting the additive model #
snp.col <- grep("rs", colnames(dat))
output <- NULL
for (i in 1:length(snp.col))
  {
  model <- glm(dat[,"affection"] ~ dat[,snp.col[i]], family="binomial")
  out.tmp1 <- summary(model)$coeff[grep("snp", rownames(summary(model)$coef)),]
  output <- rbind(output, out.tmp1)
  }
rownames(output) <- colnames(dat)[snp.col]


## Section 2 ##
## --------- ##
ped <- read.table("D:/GGCP_April2010/Practicals/Day3_Session4/gwas.ped", stringsAsFactors=F)
dim(ped)
ped[1:10,1:10]

map <- read.table("D:/GGCP_April2010/Practicals/Day3_Session4/gwas.map", stringsAsFactors=F)
dim(map)
head(map)

# 2g. Manhattan plotting #
dat.plot <- read.table("D:/GGCP_April2010/Practicals/Day3_Session4/gwas_diabetes_adj_gender.assoc.logistic", header=T,stringsAsFactors=F)
dat.snp.plot <- dat.plot[grep("ADD", dat.plot[,"TEST"]),]
plot(dat.snp.plot$BP/10^6, -log10(dat.snp.plot$P), main="Diabetes, gender adjusted",
     xlab="base position (in Mb)", ylab="-log10 p-values", pch=19)
     
# 2h. QQ-plot/PP-plot #
dat.plot <- read.table("D:/GGCP_April2010/Practicals/Day3_Session4/gwas_diabetes_adj_gender.assoc.logistic", header=T,stringsAsFactors=F)
dat.snp.plot <- dat.plot[grep("ADD", dat.plot[,"TEST"]),]
size <- nrow(dat.snp.plot)

# QQ-plot #
chisq.exp <- qchisq(seq(from=1, to=size, by=1)/(size+1), df=1)
chisq.obs <- sort((dat.snp.plot$STAT)^2)
plot(chisq.exp, chisq.obs, xlab="expected chisq statistics", pch=20,
     ylab="observed chisq statistics", main="Diabetes, gender adjusted")
abline(a=0, b=1)
x11()

# PP-plot #
p.exp <- (-log10(qunif(seq(from=1, to=size, by=1)/(size+1))))
p.obs <- (-log10(sort(dat.snp.plot[,"P"])))
plot(p.exp, p.obs, xlab="-log(expected P values)", pch=20,
     ylab="-log(observed Pvalues)", main="Diabetes, gender adjusted")
abline(a=0, b=1)







