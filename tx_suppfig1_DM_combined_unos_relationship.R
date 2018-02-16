setwd("/Users/chenj220/Documents/ucsf/transplantation/a_donor/unos/rship")

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
idata <- read.delim("unos-LIVING_DONOR_DATA_DM_race_organ_relation-kidney-us-128407.rshipmod.txt", header = T, sep = "\t")

## split at age 25
idata.lt25 <- idata[idata$AGE_DON <= 25,]
idata.mt25 <- idata[idata$AGE_DON > 25,]

##-----------------------------------------------------------
## relationships
## BAR PLOT
x11(type="cairo")
age.rship <- ddply(.data = idata, .variables = c("AGE_DON", "LIV_DON_TY_MERGED"), .fun=summarise,
                   count <- length(AGE_DON)) 
colnames(age.rship) <- c("AGE_DON", "LIV_DON_TY_MERGED", "count")
age.rship <- age.rship %>% filter(!is.na(LIV_DON_TY_MERGED))

pmain <- ggplot(age.rship, aes(AGE_DON, count, fill = LIV_DON_TY_MERGED), xlim = c(0,80))
phisto <- geom_bar(stat = "identity")
plabels <- labs(x="Age",y="Count")
pticks <- scale_x_continuous(breaks=seq(min(age.rship$AGE_DON),max(age.rship$AGE_DON),by=1))
paxes <-  theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
ptitle <- theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend <- theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
pline <- geom_vline(xintercept=25, color="red")
## set1 usually has only 8 colors
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
## how many colors do we need
colorCount <- length(unique(age.rship$LIV_DON_TY_MERGED))
## "manually" add colors instead of scale_fill_brewer fixed brwer colors
# newcolors <- c("#E41A1C","#66628D","#66628D","#419486","#5A9D5A","#66628D","#999999","#999999","black","#999999","#999999","#999999")
pmain + phisto + plabels + pticks + paxes + ptitle + plegend + scale_fill_manual(values = getPalette(colorCount)) + pline
# pmain + phisto + plabels + pticks + paxes + ptitle + plegend + scale_fill_manual(values = newcolors) + pline
# ggsave("age_reciprship_bar_newcolors.pdf", device = "pdf")

#
# LINE
pmain <- ggplot(age.rship, aes(AGE_DON, count, color = LIV_DON_TY_MERGED), xlim = c(0,80))
phisto <- geom_line(size=1.5)
pmain + phisto + plabels + pticks + paxes + ptitle + plegend + scale_color_manual(values = getPalette(colorCount)) + pline
# pmain + phisto + plabels + pticks + paxes + ptitle + plegend + scale_color_manual(values = newcolors) + pline
# ggsave("age_reciprship_line_newcolors.pdf", device = "pdf")

#
# split by gender
## BAR
age.rship.gender <- ddply(.data = idata, .variables = c("AGE_DON", "GENDER", "LIV_DON_TY_MERGED"), .fun=summarise, count = length(AGE_DON))
age.rship.gender <- age.rship.gender %>% filter(!is.na(LIV_DON_TY_MERGED))

pmain <- ggplot(age.rship.gender, aes(AGE_DON, count, fill = LIV_DON_TY_MERGED))
phisto <- geom_bar(stat = "identity")
pfacet <- facet_grid(GENDER ~., scale = "free_y")
pmain + phisto + plabels + pticks + paxes + ptitle + plegend + pline + pfacet 
# pmain + phisto + plabels + pticks + paxes + ptitle + plegend + scale_fill_manual(values = newcolors) + pline + pfacet
# ggsave("age_gender_split_reciprship_bar.pdf", device = "pdf")

## LINE
pmain <- ggplot(age.rship.gender, aes(AGE_DON, count, color = LIV_DON_TY_MERGED), xlim = c(0,80), ylim = c(0,70))
phisto <- geom_line(size=1.5)
pfacet <- facet_grid(GENDER ~., scale = "free_y")
pmain + phisto + plabels + pticks + paxes + ptitle + plegend + scale_color_manual(values = getPalette(colorCount)) + pline + pfacet
# pmain + phisto + plabels + pticks + paxes + ptitle + plegend + scale_color_manual(values = newcolors) + pline + pfacet
# ggsave("/Users/chenj220/Documents/ucsf/transplantation/manuscript/figures/suppfig2-relationship/rship_unos_gender_suppfig2.pdf", device = "pdf", width = 28, height = 16)

#---------------------------
# boxplots
cbage <- 25

acompare <- idata %>% 
  select(AGE_DON, GENDER, LIV_DON_TY_MERGED) %>%
  group_by(LIV_DON_TY_MERGED, GENDER, AGE_DON) %>%
  summarize(count = n()) %>%
  mutate(ifcbage = ifelse(AGE_DON <= cbage, 
                       paste0("lt",cbage), 
                       paste0("mt",cbage)))
  
x11(type="cairo")

ggplot(data = acompare, aes(x=ifcbage,y=count)) +
  geom_boxplot(aes(fill=GENDER)) +
  facet_wrap(~ LIV_DON_TY_MERGED,
             scales = "free") +
  scale_fill_manual(values=c("#e41a1c", "#377eb8"))

# ggsave("/Users/chenj220/Documents/ucsf/transplantation/manuscript/figures/suppfig2-relationship/rship_unos_gender_box_suppfig2.pdf", device = "pdf")

#---------------------------
# mann whitneys test
# rship = relationship
# cbage = childbearing age
# data  = acompare; contains:
# GENDER,LIV_DON_TY_MERGED,count
mwt <- function(rship, cbage, data, bonn_n){
  lt_cbage_f <- data %>%
    filter(ifcbage == paste0("lt",cbage) & GENDER == "F") %>%
    filter(LIV_DON_TY_MERGED == rship)
  
  lt_cbage_m <- data %>%
    filter(ifcbage == paste0("lt",cbage) & GENDER == "M") %>%
    filter(LIV_DON_TY_MERGED == rship)
  
  mt_cbage_f <- data %>%
    filter(ifcbage == paste0("mt",cbage) & GENDER == "F") %>%
    filter(LIV_DON_TY_MERGED == rship)
  
  mt_cbage_m <- data %>%
    filter(ifcbage == paste0("mt",cbage) & GENDER == "M") %>%
    filter(LIV_DON_TY_MERGED == rship)
  
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
  sapply(na.omit(as.character(unique(acompare$LIV_DON_TY_MERGED))),
            function(x) mwt(x, cbage, acompare, bonn_n = 6))
  )
)

write.table(a, "/Users/chenj220/Documents/ucsf/transplantation/manuscript/figures/suppfig2-relationship/rship_unos_gender_mannwhitney_suppfig2.txt", 
            quote = F, row.names = F, sep = "\t")

mwt("Biological_BloodRelatedChild", cbage, acompare, bonn_n = 6)


##!~!~!~~~~~~~~!~!~!~!~!~!~!~!~!~!
## compare ALL male and female proportions before and after cbage


# ggplot(data = acompare, aes(x=PHI_GENDER,y=count)) +
#   geom_boxplot(aes(fill=ifcbage)) +
#   facet_wrap(~ DEMO_DNR_RELATE_COMB, 
#              scales = "free")
mwt_thresh_all <- function(cbage, data, plot=0, bon_n){
  acompare_all <- data %>% 
    select(AGE_DON, GENDER) %>%
    group_by(GENDER, AGE_DON) %>%
    summarize(count = n()) %>%
    mutate(ifcbage = ifelse(AGE_DON <= cbage, 
                            paste0("lt",cbage), 
                            paste0("mt",cbage)))
  
  if(plot){
    x11(type="cairo")
    h <- ggplot(data = acompare_all, aes(x=ifcbage,y=count)) +
      geom_boxplot(aes(fill=GENDER)) +
      scale_fill_manual(values=c("#e41a1c", "#377eb8"))
    print(h)
    ggsave("/Users/chenj220/Documents/ucsf/transplantation/manuscript/figures/suppfig1-cutoff25/rship_unos_gender_box_suppfig1_all.pdf", device = "pdf")
  }

  
  lt_cbage_f <- acompare_all %>%
    filter(ifcbage == paste0("lt",cbage) & GENDER == "F")
  
  lt_cbage_m <- acompare_all %>%
    filter(ifcbage == paste0("lt",cbage) & GENDER == "M")
  
  mt_cbage_f <- acompare_all %>%
    filter(ifcbage == paste0("mt",cbage) & GENDER == "F") 
  
  mt_cbage_m <- acompare_all %>%
    filter(ifcbage == paste0("mt",cbage) & GENDER == "M") 
  
  df <- data.frame(
    cb_age = cbage,
    lt_p = wilcox.test(lt_cbage_f$count,lt_cbage_m$count)$p.value,
    mt_p = wilcox.test(mt_cbage_f$count,mt_cbage_m$count)$p.value,
    lt_p_corrected = min(wilcox.test(lt_cbage_f$count,lt_cbage_m$count)$p.value * bon_n,1),
    mt_p_corrected = min(wilcox.test(mt_cbage_f$count,mt_cbage_m$count)$p.value * bon_n,1)
  )
  
  return(df)
}

aa <- t(sapply(25:35, function(x) mwt_thresh_all(x, idata, bon_n=2)))

write.table(aa, "/Users/chenj220/Documents/ucsf/transplantation/manuscript/figures/suppfig1-cutoff25/rship_unos_gender_mannwhitney_suppfig1.txt", 
            quote = F, row.names = F, sep = "\t")

sapply(25,function(x) mwt_thresh_all(x, idata, 1, bon_n=2))