setwd("/Users/chenj220/Documents/ucsf/transplantation/a_donor/immport/finals.iTx/")

## my own library
source("~/R_codes/jmRlib.R")

library(RImmPort)
library(DBI)
library(sqldf)
library(plyr)
library(dplyr)
library(tidyr)
library(RMySQL)
library(ggplot2)
library(gridExtra)
library(gapminder)
library(reshape2)
library(RColorBrewer) 

##-----------------------------------------------------------
## input data DEMO SOCIAL
dm <- read.delim("final.iTx.DM.v15.txt", header = T, sep = "\t")
# %>% 
  # select(Sub_Org_Accession, AGE_TRANSPLANT, PHI_GENDER, ETHNI)
rp <- read.delim("final.iTx.recipient.v15.txt", header = T, sep = "\t")
idata <- merge(dm, rp, by = "Sub_Org_Accession")
idata <- idata %>% 
            mutate(AGE_TRANSPLANT = round(AGE_TRANSPLANT)) %>%
            mutate(DEMO_DNR_RELATE_COMB = 
  ifelse(DEMO_DNR_RELATE == "Child", "Biological,_Child", 
  ifelse(DEMO_DNR_RELATE == "Biological,_Half_Sibling" | DEMO_DNR_RELATE == "Sister_(full_or_half)" | DEMO_DNR_RELATE == "Brother_(full_or_half)" | DEMO_DNR_RELATE == "Identical_Twin", "Biological,_Full_Sibling",
  ifelse(DEMO_DNR_RELATE == "Nephew" | DEMO_DNR_RELATE == "Niece" | DEMO_DNR_RELATE == "Aunt" | DEMO_DNR_RELATE == "Uncle" | DEMO_DNR_RELATE == "Grandfather" | DEMO_DNR_RELATE == "Grandmother" | DEMO_DNR_RELATE == "Cousin", "Biological,_Other_Relative_(Specify)",
  ifelse(DEMO_DNR_RELATE == "Father" | DEMO_DNR_RELATE == "Mother", "Biological,_Parent", 
  ifelse(DEMO_DNR_RELATE == "Unknown", NA, 
  ifelse(DEMO_DNR_RELATE == "Non-Biological,_Friend" | DEMO_DNR_RELATE == "Non-Biological,_Unrelated,_Directed_Dona" | DEMO_DNR_RELATE == "Non-Biological,_Unrelated,_Non-Directed" | DEMO_DNR_RELATE == "Non-Biological,_Unrelated_Paired_Exchang" | DEMO_DNR_RELATE == "Non-Biological,_Other_Unrelated_Directed", "Non-Biological,_Unrelated", as.character(DEMO_DNR_RELATE)))))))) %>%
  filter(DEMO_DNR_RELATE != "Living_Related" & DEMO_DNR_RELATE != "Living_related" & DEMO_DNR_RELATE != "Living_unrelated" & DEMO_DNR_RELATE != "Living_Unrelated")


## split at age 25
idata.lt25 <- idata[idata$AGE_TRANSPLANT <= 25,]
idata.mt25 <- idata[idata$AGE_TRANSPLANT > 25,]

##-----------------------------------------------------------
## relationships
## BAR PLOT
x11(type="cairo")
age.rship <- ddply(.data = idata, .variables = c("AGE_TRANSPLANT", "DEMO_DNR_RELATE_COMB"), .fun=summarise,
                   count <- length(AGE_TRANSPLANT)) 
colnames(age.rship) <- c("AGE_TRANSPLANT", "DEMO_DNR_RELATE_COMB", "count")
age.rship <- age.rship %>% filter(!is.na(DEMO_DNR_RELATE_COMB))

pmain <- ggplot(age.rship, aes(AGE_TRANSPLANT, count, fill = DEMO_DNR_RELATE_COMB), xlim = c(0,80))
phisto <- geom_bar(stat = "identity")
plabels <- labs(x="Age",y="Count")
pticks <- scale_x_continuous(breaks=seq(min(age.rship$AGE_TRANSPLANT),max(age.rship$AGE_TRANSPLANT),by=1))
paxes <-  theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
ptitle <- theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend <- theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
pline <- geom_vline(xintercept=25, color="red")
## set1 usually has only 8 colors
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
## how many colors do we need
colorCount <- length(unique(age.rship$DEMO_DNR_RELATE_COMB))
## "manually" add colors instead of scale_fill_brewer fixed brwer colors
# newcolors <- c("#E41A1C","#66628D","#66628D","#419486","#5A9D5A","#66628D","#999999","#999999","black","#999999","#999999","#999999")
pmain + phisto + plabels + pticks + paxes + ptitle + plegend + scale_fill_manual(values = getPalette(colorCount)) + pline
# pmain + phisto + plabels + pticks + paxes + ptitle + plegend + scale_fill_manual(values = newcolors) + pline
# ggsave("age_reciprship_bar_newcolors.pdf", device = "pdf")

#
# LINE
pmain <- ggplot(age.rship, aes(AGE_TRANSPLANT, count, color = DEMO_DNR_RELATE_COMB), xlim = c(0,80))
phisto <- geom_line(size=1.5)
pmain + phisto + plabels + pticks + paxes + ptitle + plegend + scale_color_manual(values = getPalette(colorCount)) + pline
# pmain + phisto + plabels + pticks + paxes + ptitle + plegend + scale_color_manual(values = newcolors) + pline
# ggsave("age_reciprship_line_newcolors.pdf", device = "pdf")

#
# split by gender
## BAR
age.rship.gender <- ddply(.data = idata, .variables = c("AGE_TRANSPLANT", "PHI_GENDER", "DEMO_DNR_RELATE_COMB"), .fun=summarise, count = length(AGE_TRANSPLANT))
colnames(age.rship.gender) <- c("AGE_TRANSPLANT", "PHI_GENDER", "DEMO_DNR_RELATE_COMB", "count")
age.rship.gender <- age.rship.gender %>% filter(!is.na(DEMO_DNR_RELATE_COMB))

pmain <- ggplot(age.rship.gender, aes(AGE_TRANSPLANT, count, fill = DEMO_DNR_RELATE_COMB))
phisto <- geom_bar(stat = "identity")
pfacet <- facet_grid(PHI_GENDER ~., scale = "free_y")
pmain + phisto + plabels + pticks + paxes + ptitle + plegend + pline + pfacet 
# pmain + phisto + plabels + pticks + paxes + ptitle + plegend + scale_fill_manual(values = newcolors) + pline + pfacet
# ggsave("age_gender_split_reciprship_bar.pdf", device = "pdf")

## LINE
x11(type="cairo")
pmain <- ggplot(age.rship.gender, aes(AGE_TRANSPLANT, count, color = DEMO_DNR_RELATE_COMB), xlim = c(0,80), ylim = c(0,70))
phisto <- geom_line(size=1.5)
pfacet <- facet_grid(PHI_GENDER ~., scale = "free_y")
pmain + phisto + plabels + pticks + paxes + ptitle + plegend + scale_color_manual(values = getPalette(colorCount)) + pline + pfacet
# pmain + phisto + plabels + pticks + paxes + ptitle + plegend + scale_color_manual(values = newcolors) + pline + pfacet
# ggsave("/Users/chenj220/Documents/ucsf/transplantation/manuscript/figures/suppfig2-relationship/rship_immport_gender_suppfig2.pdf", device = "pdf", width = 28, height = 16)

#---------------------------
# boxplots
# base counts = nonbiol_friend + nonbiol_unrelated
## debug
# nrow(idata.lt25 %>%
#   filter(DEMO_DNR_RELATE_COMB == "Non-Biological,_Unrelated" |
#            DEMO_DNR_RELATE_COMB == "Non-Biological,_Friend") %>%
#   filter(PHI_GENDER == "Male")) #48
# nrow(idata.lt25 %>%
#   filter(DEMO_DNR_RELATE_COMB == "Non-Biological,_Unrelated" |
#            DEMO_DNR_RELATE_COMB == "Non-Biological,_Friend") %>%
#   filter(PHI_GENDER == "Female")) #56
# 
# nrow(nonbiol_unr_mt25_m <- idata.mt25 %>%
#   filter(DEMO_DNR_RELATE_COMB == "Non-Biological,_Unrelated" |
#            DEMO_DNR_RELATE_COMB == "Non-Biological,_Friend") %>%
#   filter(PHI_GENDER == "Male")) #565
# nrow(nonbiol_unr_mt25_f <- idata.mt25 %>%
#   filter(DEMO_DNR_RELATE_COMB == "Non-Biological,_Unrelated" |
#            DEMO_DNR_RELATE_COMB == "Non-Biological,_Friend") %>%
#   filter(PHI_GENDER == "Female")) #713

cbage <- 25

acompare <- idata %>% 
  select(AGE_TRANSPLANT, PHI_GENDER, DEMO_DNR_RELATE_COMB) %>%
  # mutate(DEMO_DNR_RELATE_COMB1 = ifelse(DEMO_DNR_RELATE_COMB == "Non-Biological,_Unrelated" | DEMO_DNR_RELATE_COMB == "Non-Biological,_Friend", "Non-Biological,_Unrelated_Friend", DEMO_DNR_RELATE_COMB)) %>%
  # select(-DEMO_DNR_RELATE_COMB) %>%
  group_by(DEMO_DNR_RELATE_COMB, PHI_GENDER, AGE_TRANSPLANT) %>%
  summarize(count = n()) %>%
  mutate(ifcbage = ifelse(AGE_TRANSPLANT <= cbage, 
                       paste0("lt",cbage), 
                       paste0("mt",cbage)))
  
x11(type="cairo")
# ggplot(data = acompare, aes(x=PHI_GENDER,y=count)) +
#   geom_boxplot(aes(fill=ifcbage)) +
#   facet_wrap(~ DEMO_DNR_RELATE_COMB, 
#              scales = "free")

ggplot(data = acompare, aes(x=ifcbage,y=count)) +
  geom_boxplot(aes(fill=PHI_GENDER)) +
  facet_wrap(~ DEMO_DNR_RELATE_COMB,
             scales = "free") +
  scale_fill_manual(values=c("#e41a1c", "#377eb8"))
# ggsave("/Users/chenj220/Documents/ucsf/transplantation/manuscript/figures/suppfig1-relationship/rship_immport_gender_box_suppfig2.pdf", device = "pdf")

#---------------------------
# mann whitneys test
# rship = relationship
# cbage = childbearing age
# data  = acompare; contains:
# PHI_GENDER,DEMO_DNR_RELATE_COMB,count
mwt <- function(rship, cbage, data, bonn_n){
  lt_cbage_f <- data %>%
    filter(ifcbage == paste0("lt",cbage) & PHI_GENDER == "Female") %>%
    filter(DEMO_DNR_RELATE_COMB == rship)
  
  lt_cbage_m <- data %>%
    filter(ifcbage == paste0("lt",cbage) & PHI_GENDER == "Male") %>%
    filter(DEMO_DNR_RELATE_COMB == rship)
  
  mt_cbage_f <- data %>%
    filter(ifcbage == paste0("mt",cbage) & PHI_GENDER == "Female") %>%
    filter(DEMO_DNR_RELATE_COMB == rship)
  
  mt_cbage_m <- data %>%
    filter(ifcbage == paste0("mt",cbage) & PHI_GENDER == "Male") %>%
    filter(DEMO_DNR_RELATE_COMB == rship)
  
  df <- data.frame(rship=rship,
                   lt_age=wilcox.test(lt_cbage_f$count,lt_cbage_m$count)$p.value,
                   mt_age=wilcox.test(mt_cbage_f$count,mt_cbage_m$count)$p.value,
                   lt_age_bonn_p=min(wilcox.test(lt_cbage_f$count,lt_cbage_m$count)$p.value * bonn_n,1),
                   mt_age_bonn_p=min(wilcox.test(mt_cbage_f$count,mt_cbage_m$count)$p.value * bonn_n,1),
                   stringsAsFactors = F)

  return(df)
}

a <- t(
  as.data.frame(
  sapply(na.omit(unique(acompare$DEMO_DNR_RELATE_COMB)),
            function(x) mwt(x, cbage, acompare, bonn_n = 6))
  )
)
write.table(a, "/Users/chenj220/Documents/ucsf/transplantation/manuscript/figures/suppfig2-relationship/rship_immport_gender_mannwhitney_suppfig2.txt", 
            quote = F, row.names = F, sep = "\t")

##!~!~!~~~~~~~~!~!~!~!~!~!~!~!~!~!
## compare ALL male and female proportions before and after cbage
x11(type="cairo")
# ggplot(data = acompare, aes(x=PHI_GENDER,y=count)) +
#   geom_boxplot(aes(fill=ifcbage)) +
#   facet_wrap(~ DEMO_DNR_RELATE_COMB, 
#              scales = "free")

mwt_thresh_all <- function(cbage, data, plot=0, bon_n){
  
  acompare_all <- data %>% 
    select(AGE_TRANSPLANT, PHI_GENDER) %>%
    group_by(PHI_GENDER, AGE_TRANSPLANT) %>%
    summarize(count = n()) %>%
    mutate(ifcbage = ifelse(AGE_TRANSPLANT <= cbage, 
                            paste0("lt",cbage), 
                            paste0("mt",cbage)))
  
  if(plot){
    h <- ggplot(data = acompare_all, aes(x=ifcbage,y=count)) +
    geom_boxplot(aes(fill=PHI_GENDER)) +
    scale_fill_manual(values=c("#e41a1c", "#377eb8"))
  print(h)
  ggsave("/Users/chenj220/Documents/ucsf/transplantation/manuscript/figures/suppfig1-cutoff25/rship_immport_gender_box_suppfig1_all.pdf", device = "pdf")
  }
  
  lt_cbage_f <- acompare_all %>%
    filter(ifcbage == paste0("lt",cbage) & PHI_GENDER == "Female")
  
  lt_cbage_m <- acompare_all %>%
    filter(ifcbage == paste0("lt",cbage) & PHI_GENDER == "Male")
  
  mt_cbage_f <- acompare_all %>%
    filter(ifcbage == paste0("mt",cbage) & PHI_GENDER == "Female") 
  
  mt_cbage_m <- acompare_all %>%
    filter(ifcbage == paste0("mt",cbage) & PHI_GENDER == "Male") 
  
  df <- data.frame(
    cb_age = cbage,
    lt_p = wilcox.test(lt_cbage_f$count,lt_cbage_m$count)$p.value,
    mt_p = wilcox.test(mt_cbage_f$count,mt_cbage_m$count)$p.value,
    lt_p_corrected = min(wilcox.test(lt_cbage_f$count,lt_cbage_m$count)$p.value * bon_n,1),
    mt_p_corrected = min(wilcox.test(mt_cbage_f$count,mt_cbage_m$count)$p.value * bon_n,1)
    )
  
  return(df)
}

aa <- t(sapply(25:45, function(x) mwt_thresh_all(x, idata, bon_n=2)))
write.table(aa, "/Users/chenj220/Documents/ucsf/transplantation/manuscript/figures/suppfig1-cutoff25/rship_immport_gender_mannwhitney_suppfig1.txt", 
            quote = F, row.names = F, sep = "\t")


sapply(25,function(x) mwt_thresh_all(x, idata, 1, bon_n=2))
