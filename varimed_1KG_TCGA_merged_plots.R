setwd("/Users/jiemingchen/Documents/varimed/pcawg/dz_risk_var_varimed_staging_LR_final_ext_sex_eth_spop_zm/merge_zm")

## my own library
source("/Users/jiemingchen/R_codes/jmRlib.R")
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(dplyr)
library(reshape2)

## input combined risks for 1KGp3
## note that LR = exp(sum(log(LR))) of all SNPs with given dz from same sample; LR_max = SNP with max abs logLR
LR.1kg = read.delim("combined_dz_risk_1000GP_spop_zm.txt", header = T, sep = "\t", stringsAsFactors = FALSE, na.strings = "") ## curr one
LR.1kg.final = LR.1kg %>% mutate(LLR = log10(LR), LLR_max = log10(LR_max))
LR.1kg.final$dataset = "1KGP3"

## input combined risks for ICGC_TCGA
LR.cancer = read.delim("combined_dz_risk_ICGC_TCGA_spop_zm_histology_m.txt", header = T, sep = "\t", stringsAsFactors = FALSE, na.strings = "")
LR.cancer.final = LR.cancer %>% mutate(LLR = log10(LR), LLR_max = log10(LR_max))
LR.cancer.final$dataset = "ICGC_TCGA"

### some stats
## number of cancer types in ICGC_TCGA patients
cancertypes = sort(unique(LR.cancer.final$histology_abbreviation_m))

## number of broad_phenotypes
dz.1kg = sort(unique(LR.1kg.final$broad_phenotype))
dz.tcga = sort(unique(LR.cancer.final$broad_phenotype))



# Preprocessing datasets and calculate Mann Whitney test p values (unadj) and BH (adj)
## Preprocessing datasets and calculate Mann Whitney test p values (unadj) and BH (adj) ####
## this function preprocesses the datasets into 3 matrices FOR ONE POPULATION
## Loop Mann-Whitney test for TCGA vs 1KGp3 (will take a while... about 6 min)
## 1: num(dz)-by-num(cancertypes) matrix of MW p values (unadj)
## 2: num(dz)-by-num(cancertypes) matrix of size of subsets (note subsets<10 of size are NA)
## 3: melted down version of 1 with 2 columns of dz and cancertypes as primary keys for p values (unadj)
## 4: melted down version subsetted by cancertypes to match dz_cancers
## input requires columns LLR and histology_abbreviation_m and broad_phenotype
## tflags 1: MW unadjusted p values
##        2: size of datasets
##        3: numSNPs (mean num over all individuals in each dataset, risk + protective)
preprocess_mat_by_pop <- function(pop, dzes, cancertypes, cancerdata, ref1kgdata, LLRcol, tflag=0)
{
  ## nested function to make an array of dz by cancer for apply
  mwp <- function(cancerdata, ref1kgdata, pop, LLRcol, cancertype, dz, tflag=0)
  {
    tcga = subset(cancerdata, eval(parse(text=pop)) & broad_phenotype == dz & histology_abbreviation_m == cancertype)
    
    onekg = subset(ref1kgdata, eval(parse(text=pop)) & broad_phenotype == dz)
    
    # if tflag == 0, unadjusted MW pvalues
    # arbitrary min datapoint of 10 in either dataset to compute MW test
    # 2.sided unadjusted
    # tflag = troubleshooting flag
    if (tflag == 0)
    {
      return(ifelse(nrow(tcga) > 10 & nrow(onekg) > 10, wilcox.test(tcga[,LLRcol], onekg[,LLRcol])$p.value, NA))
    }
    else if (tflag == 1)
    {
      ## tflag == 1, number of individuals in each dataset
      return(paste("1KGP3:", nrow(onekg),"|TCGA:", nrow(tcga)))
    }
    else
    {
      ## tflag == 2, mean over individuals of number of SNPs = numRiskAlleles + numProtectiveAlleles
      numSNPs.tcga = mean(tcga[,"SNP_risk"] + tcga[,"SNP_protective"])
      numSNPs.1kgp3 = mean(onekg[,"SNP_risk"] + onekg[,"SNP_protective"])
      
      return(paste("1KGP3:", numSNPs.1kgp3,"|TCGA:", numSNPs.tcga))
    }
    
  }
  
  # produce 2 cancertypes-by-dz matrices: 
  # if tflag == 0, unadjusted MW pvalues
  # if tflag == 1, size of datasets (for troubleshooting) 
  if(tflag == 0)
  {
    mat.pval = as.data.frame(sapply(cancertypes, function(i) sapply(dzes, function(j) mwp(cancerdata, ref1kgdata, pop, "LLR", i, j, tflag))))
    
    ## (3) melt pval matrix for heatmap plotting, by histology, 
    ## + unadj pval + BH-adj p val
    mat.pval2 = cbind(rownames(mat.pval), mat.pval)
    colnames(mat.pval2)[1] = "broad_phenotype"
    mat.pval.m = melt(mat.pval2, variable.name = "histology_abbreviation_m", value.name = "LLR.p", id.vars = "broad_phenotype")
    
    ## BH-adj p values
    mat.pval.m$LLR.p.adj = p.adjust(mat.pval.m$LLR.p, method = "BH") ## n = 2240 (excluding NAs)
    
    mat.pval.m$rank = rank(mat.pval.m$LLR.p)
    
    ## subset1: cancer match subset
    mat.ss1.cancer.match = mat.pval.m %>% subset(broad_phenotype == "Breast_cancer" |
                                                   broad_phenotype == "Colorectal_cancer" |
                                                   broad_phenotype == "Esophageal_cancer" |
                                                   broad_phenotype == "Renal_cell_cancer" |
                                                   broad_phenotype == "Renal_cell_carcinoma" |
                                                   broad_phenotype == "HCV-induced_hepatocellular_carcinoma" |
                                                   broad_phenotype == "HBV-induced_hepatocellular_carcinoma" |
                                                   broad_phenotype == "Hepatocellular_carcinoma" |
                                                   broad_phenotype == "Lung_adenocarcinoma" |
                                                   broad_phenotype == "Lung_cancer" |
                                                   broad_phenotype == "Non-Small_cell_lung_cancer" |
                                                   broad_phenotype == "Lung_cancer" |
                                                   broad_phenotype == "Squamous_cell_carcinoma_of_lungs" |
                                                   broad_phenotype == "Follicular_lymphoma" |
                                                   broad_phenotype == "Chronic_lymphocytic_leukemia" |
                                                   broad_phenotype == "Myeloproliferative_disorders" |
                                                   broad_phenotype == "Ovarian_cancer" |
                                                   broad_phenotype == "Pancreatic_cancer" |
                                                   broad_phenotype == "Prostate_cancer" |
                                                   broad_phenotype == "Melanoma" |
                                                   broad_phenotype == "Gastric_cancer" |
                                                   broad_phenotype == "Papillary_thyroid_cancer" |
                                                   broad_phenotype == "Thyroid_cancer")
    ## return
    return(list(mat.pval, mat.pval.m, mat.ss1.cancer.match))
  }
  else if(tflag == 1)
  {
    ## tflag == 1, number of individuals in each dataset
    mat.nums = as.data.frame(sapply(cancertypes, function(i) sapply(dzes, function(j) mwp(cancerdata, ref1kgdata, pop, "LLR", i, j, tflag)))) ## debug
    
    return(mat.nums)
  }
  else 
  {
    ## tflag == 2, mean over individuals of number of SNPs = numRiskAlleles + numProtectiveAlleles
    mat.numSNPs = as.data.frame(sapply(cancertypes, function(i) sapply(dzes, function(j) mwp(cancerdata, ref1kgdata, pop, "c(SNP_risk, SNP_protective)", i, j, tflag)))) 
    
    return(mat.numSNPs)
  }
}


# ## EUR only ~~~~~~~~~~~~~~~~~~~~~~~~
# ## loop
# ## user  system elapsed
# ## 321.887  40.770 364.277
# system.time({
#   EUR.procdata = preprocess_mat_by_pop("population == \"EUR\"", dz.tcga, cancertypes, LR.cancer.final, LR.1kg.final, "LLR", tflag=0)
# })
# EUR.mat.pval = as.data.frame(EUR.procdata[1])
# EUR.mat.pval.m = as.data.frame(EUR.procdata[2])
# EUR.cancer.match.ss1 = as.data.frame(EUR.procdata[3])
# EUR.mat.nums = preprocess_mat_by_pop("population == \"EUR\"", dz.tcga, cancertypes, LR.cancer.final, LR.1kg.final, "LLR", tflag=1) ## debug
# EUR.mat.numSNPs = preprocess_mat_by_pop("population == \"EUR\"", dz.tcga, cancertypes, LR.cancer.final, LR.1kg.final, "c(SNP_risk, SNP_protective)", tflag=2)
# 
# ## EAS only ~~~~~~~~~~~~~~~~~~~~~~~~
# ## loop
# ## user  system elapsed
# ## 321.887  40.770 364.277
# system.time({
#   EAS.procdata = preprocess_mat_by_pop("population == \"EAS\"", dz.tcga, cancertypes, LR.cancer.final, LR.1kg.final, "LLR", tflag=0)
# })
# EAS.mat.pval = as.data.frame(EAS.procdata[1])
# EAS.mat.pval.m = as.data.frame(EAS.procdata[2])
# EAS.cancer.match.ss1 = as.data.frame(EAS.procdata[3])
# EAS.mat.nums = preprocess_mat_by_pop("population == \"EAS\"", dz.tcga, cancertypes, LR.cancer.final, LR.1kg.final, "LLR", tflag=1) ## debug
# EAS.mat.numSNPs = preprocess_mat_by_pop("population == \"EAS\"", dz.tcga, cancertypes, LR.cancer.final, LR.1kg.final, "c(SNP_risk, SNP_protective)", tflag=2)
# 
# ## AMR only ~~~~~~~~~~~~~~~~~~~~~~~~
# ## loop
# ## user  system elapsed
# ## 321.887  40.770 364.277
# system.time({
#   AMR.procdata = preprocess_mat_by_pop("population == \"AMR\"", dz.tcga, cancertypes, LR.cancer.final, LR.1kg.final, "LLR", tflag=0)
# })
# AMR.mat.pval = as.data.frame(AMR.procdata[1])
# AMR.mat.pval.m = as.data.frame(AMR.procdata[2])
# AMR.cancer.match.ss1 = as.data.frame(AMR.procdata[3])
# AMR.mat.nums = preprocess_mat_by_pop("population == \"AMR\"", dz.tcga, cancertypes, LR.cancer.final, LR.1kg.final, "LLR", tflag=1) ## debug
# AMR.mat.numSNPs = preprocess_mat_by_pop("population == \"AMR\"", dz.tcga, cancertypes, LR.cancer.final, LR.1kg.final, "c(SNP_risk, SNP_protective)", tflag=2)
# 
# ## AFR only ~~~~~~~~~~~~~~~~~~~~~~~~
# ## loop
# ## user  system elapsed
# ## 321.887  40.770 364.277
# system.time({
#   AFR.procdata = preprocess_mat_by_pop("population == \"AFR\"", dz.tcga, cancertypes, LR.cancer.final, LR.1kg.final, "LLR", tflag=0)
# })
# AFR.mat.pval = as.data.frame(AFR.procdata[1])
# AFR.mat.pval.m = as.data.frame(AFR.procdata[2])
# AFR.cancer.match.ss1 = as.data.frame(AFR.procdata[3])
# AFR.mat.nums = preprocess_mat_by_pop("population == \"AFR\"", dz.tcga, cancertypes, LR.cancer.final, LR.1kg.final, "LLR", tflag=1) ## debug
# AFR.mat.numSNPs = preprocess_mat_by_pop("population == \"AFR\"", dz.tcga, cancertypes, LR.cancer.final, LR.1kg.final, "c(SNP_risk, SNP_protective)", tflag=2)
# 
# ## SAS only ~~~~~~~~~~~~~~~~~~~~~~~~
# ## loop
# ## user  system elapsed
# ## 321.887  40.770 364.277
# system.time({
#   SAS.procdata = preprocess_mat_by_pop("population == \"SAS\"", dz.tcga, cancertypes, LR.cancer.final, LR.1kg.final, "LLR", tflag=0)
# })
# SAS.mat.pval = as.data.frame(SAS.procdata[1])
# SAS.mat.pval.m = as.data.frame(SAS.procdata[2])
# SAS.cancer.match.ss1 = as.data.frame(SAS.procdata[3])
# SAS.mat.nums = preprocess_mat_by_pop("population == \"SAS\"", dz.tcga, cancertypes, LR.cancer.final, LR.1kg.final, "LLR", tflag=1) ## debug
# SAS.mat.numSNPs = preprocess_mat_by_pop("population == \"SAS\"", dz.tcga, cancertypes, LR.cancer.final, LR.1kg.final, "c(SNP_risk, SNP_protective)", tflag=2)


# Plot matrices
#################################################################
# Heatmap of p.values (adj) of cancertypes vs diseases in VariMed
## ggplot cant do heatmap v well
## ggplot cant border cells v well

library(data.table)
plotmatrix <- function(mat, colp, xfontsize, yfontsize, colors, pop)
{
  
  
  ## discretize/categorize p values into 3 categories with new column p.cat
  mat$p.cat = ifelse(mat[,colp] <= 0.01, "0-0.01", 
                     ifelse(mat[,colp] > 0.01 & mat[,colp] <=0.05, "0.01-0.05", 
                            ifelse(mat[,colp] > 0.05 & mat[,colp] <=0.1, "0.05-0.1", "0.1-")))
  
  ## change >0.1 to NA
  # mat[,colp] = ifelse(mat[,colp] > 0.1, NA, mat[,colp])
  # mat = mat[complete.cases(mat),]
  
  
  
  ## heatmap 1 -- everything
  pl = ggplot(mat, aes(histology_abbreviation_m, broad_phenotype)) + 
    geom_tile(aes(fill = p.cat), colour = "white") +
    theme(legend.position = "none", 
          axis.text.x = element_text(size = xfontsize, angle = 330, 
                                     hjust = 0, color = "black"),
          axis.text.y = element_text(size = yfontsize)) + 
    theme(legend.position="right") + 
    scale_fill_manual(values = colors, na.value = "white") +
    labs(y=paste(pop,"_broad_phenotype"),x="histology_abbreviation_m")
  # +
  # theme(panel.border=element_rect(fill = NA, colour=alpha('black', 0.5), size=5)) 
  
  pl
}

## EUR
x11(type="cairo")
plotmatrix(EUR.mat.pval.m, colp="LLR.p.adj", xfontsize=12, yfontsize=8, colors=c("black","#3794bf","#df8640","gray90"), "EUR")
# ggsave("m_EUR_alldz_allcancertypes_zm.pdf", device = "pdf")
x11(type="cairo")
plotmatrix(EUR.cancer.match.ss1, colp="LLR.p.adj", xfontsize=18, yfontsize=18, colors=c("black","gray90","#3794bf","#df8640"), "EUR")
# ggsave("m_EUR_matched_allcancertypes_zm.pdf", device = "pdf")

## plot only complete and non-NA cases
x11(type="cairo")
EUR.mat.pval.m.new = EUR.mat.pval.m
EUR.mat.pval.m.new[,"LLR.p.adj"][EUR.mat.pval.m.new$LLR.p.adj>0.1] = NA
plotmatrix(EUR.mat.pval.m.new[complete.cases(EUR.mat.pval.m.new),], colp="LLR.p.adj", xfontsize=18, yfontsize=18, colors=c("black","#3794bf","#df8640","gray90"), "EUR")
# ggsave("m_EUR_complete_cases_zm.pdf", device = "pdf")

####### bubble plot for complete cases EUR  #######
cc.eur = EUR.mat.pval.m.new[complete.cases(EUR.mat.pval.m.new),]
cc.eur$p.cat = ifelse(cc.eur[,"LLR.p.adj"] <= 0.01, "0-0.01", 
                      ifelse(cc.eur[,"LLR.p.adj"] > 0.01 & cc.eur[,"LLR.p.adj"] <=0.05, "0.01-0.05", 
                             ifelse(cc.eur[,"LLR.p.adj"] > 0.05 & cc.eur[,"LLR.p.adj"] <=0.1, "0.05-0.1", "0.1-")))
colors=c("black","#3794bf","#df8640","gray90")

cc.eur$numSNPs = unlist(mapply(function(x,y) return(data.frame(EUR.mat.numSNPs[x,y])), cc.eur$broad_phenotype, cc.eur$histology_abbreviation))
cc.eur$numSNPs = as.character(cc.eur$numSNPs)
cc.eur$onekgSNPs = 
  as.numeric(cc.eur$numSNPs %>% { gsub("1KGP3\\: ", "", .) } %>% { gsub(" \\|TCGA: .*", "", .) })
cc.eur$tcgaSNPs = 
  as.numeric(cc.eur$numSNPs %>% { gsub(".*TCGA: ", "", .) })
cc.eur$numSNPsMean = (cc.eur$onekgSNPs + cc.eur$tcgaSNPs) / 2
cc.eur$p.cat = factor(cc.eur$p.cat)


# order of cancer type [MANUAL]
cc.eur$histology_abbreviation_m = factor(cc.eur$histology_abbreviation_m, 
                                       levels=c("Prost-AdenoCA",
                                                                    "Skin-Melanoma",
                                                                    "CNS_cancer",
                                                                    "Panc_cancer",
                                                                    "Eso-AdenoCA",
                                                                    "Eso-AdenoCA",
                                                                    "Stomach-AdenoCA",
                                                                    "ColoRect-AdenoCA"))

cc.eur$broad_phenotype = factor(cc.eur$broad_phenotype,levels=rev(c("Prostate_cancer",
                                                                    "Melanoma",
                                                                    "Glioma",
                                                                    "Biliary_liver_cirrhosis",
                                                                    "Narcolepsy",
                                                                    "Wilms'_tumor",
                                                                    "Hair_color",
                                                                    "Melanoma",
                                                                    "Behcet's_disease",
                                                                    "Hair_color")))


write.table(cc.eur, "cc.eur.txt", sep="\t", quote=F)

## p values for more, less and 2 sided
mwp <- function(cancerdata, ref1kgdata, pop, LLRcol, cancertype, dz)
{
  tcga = subset(cancerdata, eval(parse(text=pop)) & broad_phenotype == dz & histology_abbreviation_m == cancertype)
  
  onekg = subset(ref1kgdata, eval(parse(text=pop)) & broad_phenotype == dz)
  
  # if tflag == 0, unadjusted MW pvalues
  # arbitrary min datapoint of 10 in either dataset to compute MW test
  # 2.sided unadjusted
  twosided = ifelse(nrow(tcga) > 10 & nrow(onekg) > 10, wilcox.test(tcga[,LLRcol], onekg[,LLRcol])$p.value, NA)
  
  greater = ifelse(nrow(tcga) > 10 & nrow(onekg) > 10, wilcox.test(tcga[,LLRcol], onekg[,LLRcol], alternative = "greater")$p.value, NA)
  
  less = ifelse(nrow(tcga) > 10 & nrow(onekg) > 10, wilcox.test(tcga[,LLRcol], onekg[,LLRcol], alternative = "less")$p.value, NA)
  
  return(data.frame(dz, cancertype, twosided, greater, less, stringsAsFactors = FALSE))
  
}

## n = 2872 for multiple testing
cc.eur.bp = c("Prostate_cancer","Melanoma","Glioma","Biliary_liver_cirrhosis","Narcolepsy","Wilms'_tumor","Hair_color","Melanoma","Behcet's_disease","Hair_color","Prostate_cancer")
cc.eur.ha = c("Prost-AdenoCA","Skin-Melanoma","CNS_cancer","CNS_cancer","Panc_cancer","Panc_cancer","Eso-AdenoCA","Eso-AdenoCA","Stomach-AdenoCA","ColoRect-AdenoCA","CNS_cancer")
mwpp.eur = as.data.frame(t(mapply(function(x,y) mwp(LR.cancer.final, LR.1kg.final, "population==\"EUR\"", "LLR", x,y), cc.eur.ha, cc.eur.bp)))
colnames(mwpp.eur) = c("broad_phenotype","histology_abbreviation","twosided","greater","less")
mwpp.eur$risk = ifelse(as.numeric(mwpp.eur$greater) < as.numeric(mwpp.eur$less), "greater", "less")
write.table(as.matrix(mwpp.eur), "cc.eur.mwp.txt", sep="\t", quote=F)


# plot
x11(type="cairo")
ggplot(cc.eur, aes(x=histology_abbreviation, y=broad_phenotype, size=numSNPsMean, fill=p.cat, color=p.cat)) + 
  geom_point(shape = 21) +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 18, angle = 330, 
                                   hjust = 0, color = "black"),
        axis.text.y = element_text(size = 18)) + 
  theme(legend.position="right") + 
  scale_fill_manual(values = colors, na.value = "white") +
  scale_color_manual(values = colors, na.value = "white") + 
  labs(y=paste("EUR","_broad_phenotype"),x="histology_abbreviation") + 
  scale_size_area(max_size = 20) +
  theme(legend.title = element_text(size=15, face="bold"), legend.text=element_text(size=15)) 

## this add label for numbers
# + geom_text(aes(label=round(numSNPsMean)), size=5, nudge_x=0.0, nudge_y=0.6)
# ggsave("m_EUR_complete_cases_bubble_zm.pdf", device = "pdf", useDingbats=FALSE, width = 10.9, height = 9.36)

# Prost-AdenoCA and Prostate_cancer
gender.1kg = subset(read.delim("1kg-sampleinfo.txt", header=T, stringsAsFactors=FALSE, sep="\t"), select=c(sampleID, Gender))
LR.1kg.final.gender = merge(LR.1kg.final, gender.1kg, by.x="sample.id", by.y="sampleID", all.x=TRUE)
LR.1kg.final.male = subset(LR.1kg.final.gender, Gender == "male")
x11(type="cairo")
p.prostate.prostate = as.data.frame(plotviolin(LR.cancer.final, LR.1kg.final, 
                                               "population==\"EUR\"", 
                                               "histology_abbreviation==\"Prost-AdenoCA\"",
                                               dz = "Prostate_cancer"))
# ggsave("violin-prostate-prostate_zm.pdf", device = "pdf", useDingbats=FALSE)
x11(type="cairo")
p.prostate.prostate.m = as.data.frame(plotviolin(LR.cancer.final, LR.1kg.final.male, 
                                                 "population==\"EUR\"", 
                                                 "histology_abbreviation==\"Prost-AdenoCA\"",
                                                 dz = "Prostate_cancer"))
names(p.prostate.prostate.m) = c("twosided","onesided_less","onesided_greater")
p.prostate.prostate.m$twosided.adj = p.prostate.prostate.m$twosided * 
  (cc.eur[cc.eur$broad_phenotype == "Prostate_cancer" & cc.eur$histology_abbreviation == "Prost-AdenoCA","LLR.p.adj"] / 
     cc.eur[cc.eur$broad_phenotype == "Prostate_cancer" & cc.eur$histology_abbreviation == "Prost-AdenoCA","LLR.p"] * 
     cc.eur[cc.eur$broad_phenotype == "Prostate_cancer" & cc.eur$histology_abbreviation == "Prost-AdenoCA","rank"])
# ggsave("m_violin-prostate-prostate-male-only_zm.pdf", device = "pdf", useDingbats=FALSE)


####### ####### ####### ####### ####### ####### ####### 
## EAS
x11(type="cairo")
plotmatrix(EAS.mat.pval.m, colp="LLR.p.adj", xfontsize=12, yfontsize=8, colors=c("black","#3794bf","#df8640","gray90"), "EAS")
# ggsave("m_EAS_alldz_allcancertypes_zm.pdf", device = "pdf")
x11(type="cairo")
plotmatrix(EAS.cancer.match.ss1, colp="LLR.p.adj", xfontsize=15, yfontsize=15, colors=c("black","#3794bf","gray90","#df8640"), "EAS")
# ggsave("m_EAS_matched_allcancertypes_zm.pdf", device = "pdf")

## plot only complete and non-NA cases
x11(type="cairo")
EAS.mat.pval.m.new = EAS.mat.pval.m
EAS.mat.pval.m.new[,"LLR.p.adj"][EAS.mat.pval.m.new$LLR.p.adj>0.1] = NA
plotmatrix(EAS.mat.pval.m.new[complete.cases(EAS.mat.pval.m.new),], colp="LLR.p.adj", xfontsize=18, yfontsize=18, colors=c("black","#3794bf","#df8640","gray90"), "EAS")
# ggsave("m_EAS_complete_cases_zm.pdf", device = "pdf")

################################################################################
####### bubble plot for complete cases EAS #######
cc.eas = EAS.mat.pval.m.new[complete.cases(EAS.mat.pval.m.new),]
cc.eas$p.cat = ifelse(cc.eas[,"LLR.p.adj"] <= 0.01, "0-0.01", 
                      ifelse(cc.eas[,"LLR.p.adj"] > 0.01 & cc.eas[,"LLR.p.adj"] <=0.05, "0.01-0.05", 
                             ifelse(cc.eas[,"LLR.p.adj"] > 0.05 & cc.eas[,"LLR.p.adj"] <=0.1, "0.05-0.1", "0.1-")))
colors=c("black","#3794bf","#df8640","gray90")

cc.eas$numSNPs = unlist(mapply(function(x,y) return(data.frame(EAS.mat.numSNPs[x,y])), cc.eas$broad_phenotype, cc.eas$histology_abbreviation))
cc.eas$numSNPs = as.character(cc.eas$numSNPs)
cc.eas$onekgSNPs = 
  as.numeric(cc.eas$numSNPs %>% { gsub("1KGP3\\: ", "", .) } %>% { gsub(" \\|TCGA: .*", "", .) })
cc.eas$tcgaSNPs = 
  as.numeric(cc.eas$numSNPs %>% { gsub(".*TCGA: ", "", .) })
cc.eas$numSNPsMean = (cc.eas$onekgSNPs + cc.eas$tcgaSNPs) / 2
cc.eas$p.cat = factor(cc.eas$p.cat)


# order of cancer type
cc.eas$broad_phenotype = factor(cc.eas$broad_phenotype, 
                                levels=rev(c("HBV-induced_hepatocellular_carcinoma","Thyroid_cancer","Esophageal_cancer","Nasopharyngeal_carcinoma","Lung_cancer","Atopic_eczema","Celiac_disease","Inflammatory_bowel_disease","Ulcerative_colitis","Coronary_artery_disease","Duodenal_ulcer","Glaucoma","Graves'_disease","Kawasaki_disease","Multiple_sclerosis","Myocardial_infarction","Narcolepsy","Polycystic_ovary_syndrome","Primary_biliary_cirrhosis","Rheumatoid_arthritis","Schizophrenia","Systemic_lupus_erythematosus","Systemic_sclerosis","Type_1_diabetes","Type_2_diabetes","Uterine_leiomyoma","Vitiligo","Intracranial_aneurysm")))

write.table(cc.eas, "cc.eas.txt", sep="\t", quote=F)

## cc.eas.bp and cc.eas.ha have to be in order and matched
cc.eas.bp = c("HBV-induced_hepatocellular_carcinoma","Thyroid_cancer","Esophageal_cancer","Nasopharyngeal_carcinoma","Lung_cancer","Atopic_eczema","Celiac_disease","Inflammatory_bowel_disease","Ulcerative_colitis","Coronary_artery_disease","Duodenal_ulcer","Glaucoma","Graves'_disease","Kawasaki_disease","Multiple_sclerosis","Myocardial_infarction","Narcolepsy","Polycystic_ovary_syndrome","Primary_biliary_cirrhosis","Rheumatoid_arthritis","Schizophrenia","Systemic_lupus_erythematosus","Systemic_sclerosis","Type_1_diabetes","Type_2_diabetes","Uterine_leiomyoma","Vitiligo","Intracranial_aneurysm")

cc.eas.ha = c("Liver-HCC","Liver-HCC","Liver-HCC","Liver-HCC","Liver-HCC","Liver-HCC","Liver-HCC","Liver-HCC","Stomach-AdenoCA","Liver-HCC","Liver-HCC","Liver-HCC","Liver-HCC","Liver-HCC","Liver-HCC","Liver-HCC","Liver-HCC","Liver-HCC","Liver-HCC","Liver-HCC","Liver-HCC","Liver-HCC","Liver-HCC","Liver-HCC","Liver-HCC","Liver-HCC","Liver-HCC","Stomach-AdenoCA")

mwpp.eas = as.data.frame(t(mapply(function(x,y) mwp(LR.cancer.final, LR.1kg.final, "population==\"EAS\"", "LLR", x,y), cc.eas.ha, cc.eas.bp)))
colnames(mwpp.eas) = c("broad_phenotype","histology_abbreviation","twosided","greater","less")
mwpp.eas$risk = ifelse(as.numeric(mwpp.eas$greater) < as.numeric(mwpp.eas$less), "greater", "less")
write.table(as.matrix(mwpp.eas), "cc.eas.mwp.txt", sep="\t", quote=F)

# plot
ggplot(cc.eas, aes(x=histology_abbreviation_m, y=broad_phenotype, size=numSNPsMean, fill=p.cat, color=p.cat)) + 
  geom_point(shape = 21) +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 18, angle = 330, 
                                   hjust = 0, color = "black"),
        axis.text.y = element_text(size = 18)) + 
  theme(legend.position="right") + 
  scale_fill_manual(values = colors, na.value = "white") +
  scale_color_manual(values = rep("white",3), na.value = "white") + 
  labs(y=paste("EAS","_broad_phenotype"),x="histology_abbreviation") + 
  scale_size_area(max_size = 20) +
  theme(legend.title = element_text(size=15, face="bold"), legend.text=element_text(size=15)) 

## this add label for numbers
# + geom_text(aes(label=round(numSNPsMean)), size=5, nudge_x=0.0, nudge_y=0.6)
# ggsave("m_EAS_complete_cases_bubble_zm.pdf", device = "pdf", useDingbats=FALSE) 
## 10.4 x 14.4 


################################################################################



# ## AFR
# x11(type="cairo")
# plotmatrix(AFR.mat.pval.m, colp="LLR.p.adj", xfontsize=8, yfontsize=5, colors=c("gray90","black","#3794bf","#df8640"), "AFR")
# plotmatrix(AFR.cancer.match.ss1, colp="LLR.p.adj", xfontsize=5, yfontsize=7, colors=c("black","#df8640","gray90", "#3794bf"), "AFR")

x11(type="cairo")
AFR.mat.pval.m.new = AFR.mat.pval.m
AFR.mat.pval.m.new[,"LLR.p.adj"][AFR.mat.pval.m.new$LLR.p.adj>0.1] = NA
plotmatrix(AFR.mat.pval.m.new[complete.cases(AFR.mat.pval.m.new),], colp="LLR.p.adj", xfontsize=18, yfontsize=18, colors=c("black","#3794bf","#df8640","gray90"), "AFR")
# 
# ## AMR
# x11(type="cairo")
# plotmatrix(AMR.mat.pval.m, colp="LLR.p.adj", xfontsize=8, yfontsize=5, colors=c("black","#3794bf","#df8640","gray90"), "AMR")
# plotmatrix(AMR.cancer.match.ss1, colp="LLR.p.adj", xfontsize=5, yfontsize=7, colors=c("black","#df8640","gray90", "#3794bf"), "AMR")
# 
# ## SAS
# plotmatrix(SAS.mat.pval.m, colp="LLR.p.adj", xfontsize=8, yfontsize=5, colors=c("gray90", "black","#3794bf","#df8640"), "SAS")
# plotmatrix(SAS.cancer.match.ss1, colp="LLR.p.adj", xfontsize=5, yfontsize=7, colors=c("black","#df8640","gray90", "#3794bf"), "SAS")



# violin plots

plotviolin <- function(cancerdata, refdata, popparse, cancertypeparse, dz)
{
  tcga = subset(LR.cancer.final, eval(parse(text=popparse)) & eval(parse(text=cancertypeparse)), 
                select=c(sample.id, population, broad_phenotype, LLR, LLR_max, dataset))
  kgp3 = subset(LR.1kg.final, eval(parse(text=popparse)), 
                select=c(sample.id, population, broad_phenotype, LLR, LLR_max, dataset))
  
  merged = rbind(tcga, kgp3)
  
  ## mann whitney test for melanoma
  mm.1kgp3 = kgp3[kgp3$broad_phenotype==dz,]
  mm.tcga = tcga[tcga$broad_phenotype==dz,]
  mm.merged = merged[merged$broad_phenotype == dz,]
  
  print(paste(popparse,"_", cancertypeparse, "_", dz, " p.val, 2-sided"))
  jm1 = wilcox.test(mm.1kgp3$LLR, mm.tcga$LLR)$p.value
  print(paste(popparse,"_", cancertypeparse, "_", dz, " p.val, less"))
  jm2 = wilcox.test(mm.1kgp3$LLR, mm.tcga$LLR, alternative = "less")$p.value ## x < y
  print(paste(popparse,"_", cancertypeparse, "_", dz, " p.val, greater"))
  jm3 = wilcox.test(mm.1kgp3$LLR, mm.tcga$LLR, alternative = "greater")$p.value ## x > y
  
  ## plot violin and boxplot for melanoma
  pd <- position_dodge(0.9)
  pmain2 = ggplot(mm.merged, aes(x=dataset, y=LLR, fill = factor(dataset)))
  phisto2 = geom_violin()
  phisto3 = geom_boxplot(width=.1, outlier.size=0, fill="grey50", position=pd) 
  phisto4 = stat_summary(fun.y=median)
  ptitle = ggtitle(paste(gsub("population==","",popparse), " ", gsub("histology_abbreviation_m==", "", cancertypeparse)))
  plabels = labs(x=dz,y="LLR distribution")
  jm4 = pmain2 + phisto2 + phisto3 + phisto4 + 
    ptitle + plabels + theme(legend.position="none") + 
    scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1") + coord_flip() 
  # + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1, colour = "black", fill = "black")
  print(jm4)
  
  return(list(jm1,jm2,jm3))
}

## item1: 2 sided unadj pvalue, 1sided x<y, 1sided x>y
# Melanoma-Melanoma
p.melanoma.melanoma = as.data.frame(plotviolin(LR.cancer.final, LR.1kg.final, 
                                               "population==\"EUR\"", 
                                               "histology_abbreviation_m==\"Skin-Melanoma\"",
                                               dz = "Melanoma"))
names(p.melanoma.melanoma) = c("2sided","1sided_less","1sided_greater")
p.melanoma.melanoma$twosided_adj.p = EUR.mat.pval.m[EUR.mat.pval.m$broad_phenotype=="Melanoma" & EUR.mat.pval.m$histology_abbreviation_m=="Skin-Melanoma",]$LLR.p.adj
p.melanoma.melanoma$size = EUR.mat.numSNPs["Melanoma","Skin-Melanoma"]
# ggsave("m_violin-melanoma.melanoma2_zm.pdf", device = "pdf")

# melanoma-obesity
p.melanoma.obesity = as.data.frame(plotviolin(LR.cancer.final, LR.1kg.final, 
                                               "population==\"EUR\"", 
                                               "histology_abbreviation_m==\"Skin-Melanoma\"",
                                               dz = "Obesity"))
names(p.melanoma.obesity) = c("2sided","1sided_less","1sided_greater")
p.melanoma.obesity$twosided_adj.p = EUR.mat.pval.m[EUR.mat.pval.m$broad_phenotype=="Obesity" & EUR.mat.pval.m$histology_abbreviation_m=="Skin-Melanoma",]$LLR.p.adj
p.melanoma.obesity$size = EUR.mat.numSNPs["Obesity","Skin-Melanoma"]


# Prost-AdenoCA and Prostate_cancer
x11(type="cairo")
p.prostate.prostate = as.data.frame(plotviolin(LR.cancer.final, LR.1kg.final, 
                                              "population==\"EUR\"", 
                                              "histology_abbreviation_m==\"Prost-AdenoCA\"",
                                              dz = "Prostate_cancer"))
# ggsave("m_violin-prostate-prostate_zm.pdf", device = "pdf")


########################
# histograms: compare (1KGp3 EUR) vs (ICGC_TCGA EUR melanoma patients) LRs for ALL VariMed diseases
tcga.EUR.melanoma = subset(LR.cancer.final, population == "EUR" & histology_abbreviation_m == "Skin-Melanoma", 
                           select=c(sample.id, population, broad_phenotype, LLR, LLR_max, dataset))
kgp3.EUR.melanoma = subset(LR.1kg.final, population == "EUR", 
                           select=c(sample.id, population, broad_phenotype, LLR, LLR_max, dataset))

merged.EUR.melanoma = rbind(tcga.EUR.melanoma, kgp3.EUR.melanoma)

## plotting
pmain = ggplot(tcga.EUR.melanoma[tcga.EUR.melanoma$broad_phenotype %in% dz.tcga[1:10],], aes(x=broad_phenotype, y=LLR))
phisto = geom_boxplot()
ptitle = ggtitle("Skin-Melanoma")
# pfacet = facet_wrap( ~ broad_phenotype, scales="free", ncol=1) 
plabels = labs(x="broad phenotype",y="LLR distribution")
# paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
# axis.title.y = element_text(face = "bold",colour = "black", size = 20),
# axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))

pmain + phisto + ptitle + plabels + scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1") + coord_flip() + geom_jitter(height = 0, width = 0.1)

## new
by=10
for (i in seq(1,length(dz.tcga),by=by))
{
  pmain2 = ggplot(merged.EUR.melanoma[merged.EUR.melanoma$broad_phenotype %in% dz.tcga[i:(i+by-1)],], 
                  aes(x=broad_phenotype, y=LLR, fill = factor(dataset)))
  phisto2 = geom_boxplot(width=0.7, outlier.shape=3) ## shape 3 = '+'
  j = pmain2 + phisto2 + ptitle + plabels + 
    scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1") + coord_flip() 
  print(j)
}



# Mann-Whitney tests & violin plots for TCGA 'Skin-Melanoma' EUR patients for LRs for 'Melanoma' (pos) and 'Obesity' (neg) and 'Renal_cell_cancer'
## mann whitney test for melanoma
mm.1kgp3 = kgp3.EUR.melanoma[kgp3.EUR.melanoma$broad_phenotype=="Melanoma",]
mm.tcga = tcga.EUR.melanoma[tcga.EUR.melanoma$broad_phenotype=="Melanoma",]
mm.merged = merged.EUR.melanoma[merged.EUR.melanoma$broad_phenotype == "Melanoma",]

print("Melanoma-Melanoma p.val, 2-sided")
wilcox.test(mm.1kgp3$LLR, mm.tcga$LLR)$p.value
print("Melanoma-Melanoma p.val, less")
wilcox.test(mm.1kgp3$LLR, mm.tcga$LLR, alternative = "less")$p.value ## x < y
print("Melanoma-Melanoma p.val, greater")
wilcox.test(mm.1kgp3$LLR, mm.tcga$LLR, alternative = "greater")$p.value ## x > y

## plot violin and boxplot for melanoma
pd <- position_dodge(0.9)
pmain2 = ggplot(mm.merged, aes(x=dataset, y=LLR, fill = factor(dataset)))
phisto2 = geom_violin()
phisto3 = geom_boxplot(width=.1, outlier.size=0, fill="grey50", position=pd) 
phisto4 = stat_summary(fun.y=median)
plabels = labs(x="Melanoma",y="LLR distribution")
pmain2 + phisto2 + phisto3 + phisto4 + 
  ptitle + plabels + theme(legend.position="none") + 
  scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1") + coord_flip()

#---------

## mann whitney test for Renal_cell_cancer
mr.1kgp3 = kgp3.EUR.melanoma[kgp3.EUR.melanoma$broad_phenotype=="Renal_cell_cancer",]
mr.tcga = tcga.EUR.melanoma[tcga.EUR.melanoma$broad_phenotype=="Renal_cell_cancer",]
mr.merged = merged.EUR.melanoma[merged.EUR.melanoma$broad_phenotype == "Renal_cell_cancer",]

print("Melanoma-renal_cell_cancer p.val, 2-sided")
wilcox.test(mr.1kgp3$LLR, mr.tcga$LLR)$p.value
print("Melanoma-renal_cell_cancer p.val, less")
wilcox.test(mr.1kgp3$LLR, mr.tcga$LLR, alternative = "less")$p.value ## x < y
print("Melanoma-renal_cell_cancer p.val, greater")
wilcox.test(mr.1kgp3$LLR, mr.tcga$LLR, alternative = "greater")$p.value ## x > y

## plot violin and boxplot for Renal_cell_cancer
pd <- position_dodge(0.9)
pmain2 = ggplot(mr.merged, aes(x=dataset, y=LLR, fill = factor(dataset)))
phisto2 = geom_violin()
phisto3 = geom_boxplot(width=.1, outlier.size=0, fill="grey50", position=pd) 
phisto4 = stat_summary(fun.y=median)
plabels = labs(x="Renal_cell_cancer",y="LLR distribution")
pmain2 + phisto2 + phisto3 + phisto4 + 
  ptitle + plabels + theme(legend.position="none") + 
  scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1") + coord_flip()

## mann whitney test for Obesity
mo.1kgp3 = kgp3.EUR.melanoma[kgp3.EUR.melanoma$broad_phenotype=="Obesity",]
mo.tcga = tcga.EUR.melanoma[tcga.EUR.melanoma$broad_phenotype=="Obesity",]
mo.merged = merged.EUR.melanoma[merged.EUR.melanoma$broad_phenotype == "Obesity",]


print("Melanoma-Obesity p.val, 2-sided")
wilcox.test(mo.1kgp3$LLR, mo.tcga$LLR)$p.value
print("Melanoma-Obesity p.val, less")
wilcox.test(mo.1kgp3$LLR, mo.tcga$LLR, alternative = "less")$p.value ## x < y
print("Melanoma-Obesity p.val, greater")
wilcox.test(mo.1kgp3$LLR, mo.tcga$LLR, alternative = "greater")$p.value ## x > y

## plot violin and boxplots for obesity
pmain2 = ggplot(mo.merged, aes(x=dataset, y=LLR, fill = factor(dataset)))
phisto2 = geom_violin()
phisto3 = geom_boxplot(width=.1, outlier.size=0, fill="grey50", position=pd) 
phisto4 = stat_summary(fun.y=median)
plabels = labs(x="Obesity",y="LLR distribution")
pmain2 + phisto2 + phisto3 + phisto4 + 
  ptitle + plabels + theme(legend.position="none") + 
  scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1") + coord_flip()




