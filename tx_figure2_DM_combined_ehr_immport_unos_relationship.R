setwd("/Users/jiemingchen/Documents/transplantation/a_donor/immport/finals.iTx")

## my own library
source("/Users/jiemingchen/R_codes/jmRlib.R")

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
dm <- read.delim("final.iTx.DM.v15.txt", header = T, sep = "\t") %>% select(Sub_Org_Accession, AGE_TRANSPLANT, PHI_GENDER)
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
  ifelse(DEMO_DNR_RELATE == "Non-Biological,_Unrelated,_Directed_Dona" | DEMO_DNR_RELATE == "Non-Biological,_Unrelated,_Non-Directed" | DEMO_DNR_RELATE == "Non-Biological,_Unrelated_Paired_Exchang" | DEMO_DNR_RELATE == "Non-Biological,_Other_Unrelated_Directed", "Non-Biological,_Unrelated", as.character(DEMO_DNR_RELATE)))))))) %>%
  filter(DEMO_DNR_RELATE != "Living_Related" & DEMO_DNR_RELATE != "Living_related" & DEMO_DNR_RELATE != "Living_unrelated" & DEMO_DNR_RELATE != "Living_Unrelated")


## split at age 25
idata.lt25 <- idata[idata$AGE_TRANSPLANT <= 25,]
idata.mt25 <- idata[idata$AGE_TRANSPLANT > 25,]

##-----------------------------------------------------------
## relationships
## BAR PLOT
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
age.rship.gender <- ddply(.data = idata, .variables = c("AGE_TRANSPLANT", "PHI_GENDER", "DEMO_DNR_RELATE_COMB"), .fun=summarise, count = length(AGE_TRANSPLANT), xlim = c(0,80))
colnames(age.rship.gender) <- c("AGE_TRANSPLANT", "PHI_GENDER", "DEMO_DNR_RELATE_COMB", "count")
age.rship.gender <- age.rship.gender %>% filter(!is.na(DEMO_DNR_RELATE_COMB))

pmain <- ggplot(age.rship.gender, aes(AGE_TRANSPLANT, count, fill = DEMO_DNR_RELATE_COMB))
phisto <- geom_bar(stat = "identity")
pfacet <- facet_grid(PHI_GENDER ~., scale = "free_y")
pmain + phisto + plabels + pticks + paxes + ptitle + plegend + pline + pfacet 
# pmain + phisto + plabels + pticks + paxes + ptitle + plegend + scale_fill_manual(values = newcolors) + pline + pfacet
# ggsave("age_gender_split_reciprship_bar.pdf", device = "pdf")

## LINE
pmain <- ggplot(age.rship.gender, aes(AGE_TRANSPLANT, count, color = DEMO_DNR_RELATE_COMB), xlim = c(0,80), ylim = c(0,70))
phisto <- geom_line(size=1.5)
pfacet <- facet_grid(~ PHI_GENDER, scale = "free_y")
pmain + phisto + plabels + pticks + paxes + ptitle + plegend + scale_color_manual(values = getPalette(colorCount)) + pline + pfacet
# pmain + phisto + plabels + pticks + paxes + ptitle + plegend + scale_color_manual(values = newcolors) + pline + pfacet
# ggsave("age_gender_split_reciprship_line_newcolors.pdf", device = "pdf")


#
# split by gender and ethnicity; no age


##-----------------------------------------------------------

## gender plots
b1.d=ddply(.data = a.d, .variables = c("AGE_TRANSPLANT","PHI_GENDER"), .fun=summarise, count = length(AGE_TRANSPLANT))

pmain = ggplot(b1.d, aes(AGE_TRANSPLANT, count, color = PHI_GENDER))
phisto = geom_line(size=1.5)
plabels = labs(x="Age",y="Count")
pticks = scale_x_continuous(breaks=seq(min(b3.d$AGE_TRANSPLANT),max(b3.d$AGE_TRANSPLANT),by=1))
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))

medM = median(idata[(idata$PHI_GENDER=="Male"),]$AGE_TRANSPLANT)
medF = median(idata[(idata$PHI_GENDER=="Female"),]$AGE_TRANSPLANT)

pmain + phisto + plabels + pticks + paxes + ptitle + plegend + scale_color_brewer(palette = "Set1")
# ggsave("gender_line_relive.pdf", device = "pdf")


## ethnicity/race; gender split by race 
b3.d=ddply(.data = idata, .variables = c("AGE_TRANSPLANT","DEMO_RACE1","PHI_GENDER"), .fun=summarise, count = length(AGE_TRANSPLANT))
total = dcast(b3.d, AGE_TRANSPLANT ~ PHI_GENDER , sum, value.var="count")
total = melt(total, id=c("AGE_TRANSPLANT"), value.name = c("count"))
total = plyr::rename(total, c("variable"="PHI_GENDER"))
total$DEMO_RACE1 = "Total"
b3.d.t = rbind(b3.d,total)

pmain = ggplot(b3.d.t, aes(AGE_TRANSPLANT, count, color = PHI_GENDER))
phisto = geom_line(size=1.5)
pfacet = facet_grid(DEMO_RACE1 ~., scale = "free_y") 
plabels = labs(x="Age",y="Count")
pticks = scale_x_continuous(breaks=seq(min(b3.d$AGE_TRANSPLANT),max(b3.d$AGE_TRANSPLANT),by=1))
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
pmain + phisto + pfacet + plabels + pticks + paxes + ptitle + plegend + scale_color_brewer(palette = "Set1")
# ggsave("gender_race_split_after_line.pdf", device = "pdf")



###- - - - -- - - - - - -- ------- - -- - --- -- ------- - -- --- -- --- --
### split 25, split gender
lt25.age.eth = ddply(.data = idata.lt25, .variables = c("AGE_TRANSPLANT", "DEMO_RACE1", "PHI_GENDER"), .fun = summarise, count = length(AGE_TRANSPLANT))
total = dcast(lt25.age.eth, AGE_TRANSPLANT ~ PHI_GENDER , sum, value.var="count")
total = melt(total, id=c("AGE_TRANSPLANT"), value.name = c("count"))
total = plyr::rename(total, c("variable"="PHI_GENDER"))
total$DEMO_RACE1 = "Total"
lt25.age.eth.t = rbind(lt25.age.eth,total)
pmain = ggplot(lt25.age.eth.t, aes(x=PHI_GENDER, y=count, fill=PHI_GENDER)) ## remove unknown sex
phisto = geom_bar(stat="identity")
pfacet = facet_grid(DEMO_RACE1 ~., scale = "free_y") 
plabels = labs(x="Gender",y="Count (lt25)")
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
h = pmain + phisto + pfacet + plabels + paxes + ptitle + plegend + scale_fill_brewer(palette = "Set1")
hh = h %+% lt25.age.eth.t

## mt25
mt25.age.eth = ddply(.data = idata.mt25, .variables = c("AGE_TRANSPLANT", "DEMO_RACE1", "PHI_GENDER"), .fun = summarise, count = length(AGE_TRANSPLANT))
total = dcast(mt25.age.eth, AGE_TRANSPLANT ~ PHI_GENDER , sum, value.var="count")
total = melt(total, id=c("AGE_TRANSPLANT"), value.name = c("count"))
total = plyr::rename(total, c("variable"="PHI_GENDER"))
total$DEMO_RACE1 = "Total"
mt25.age.eth.t = rbind(mt25.age.eth, total)
plabels = labs(x="Gender",y="Count (mt25)")
i = pmain + phisto + pfacet + plabels + paxes + ptitle + plegend + scale_fill_brewer(palette = "Set1")
ii = i %+% mt25.age.eth.t

grid.arrange(hh,ii, ncol=2)
g <- arrangeGrob(hh, ii, ncol=2) #generates g
# ggsave("gender_race_split_after_line_lt25_mt25.pdf", g, device = "pdf")


###- - - - -- - - - - - -- ------- - -- - --- -- ------- - -- --- -- --- --
### split 25, split gender, reciprship 
lt25.sex.rship.race = ddply(.data = idata.lt25, .variables = c("PHI_GENDER", "DEMO_DNR_RELATE", "DEMO_RACE1"), .fun = summarise, count = length(PHI_GENDER))
total2 = dcast(lt25.sex.rship.race, PHI_GENDER ~ DEMO_DNR_RELATE, sum, value.var="count")
total2 = melt(total2, id=c("PHI_GENDER"), value.name = c("count"))
total2 = plyr::rename(total2, c("variable"="DEMO_DNR_RELATE"))
total2$DEMO_RACE1 = "Total"
lt25.sex.rship.race.t = rbind(lt25.sex.rship.race,total2)
pmain = ggplot(lt25.sex.rship.race.t, aes(x=PHI_GENDER, y=count, fill=DEMO_DNR_RELATE)) ## remove unknown sex
phisto = geom_bar(stat="identity")
pfacet = facet_grid(DEMO_RACE1 ~., scale = "free_y") 
plabels = labs(x="Gender",y="Count (lt25)")
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
h = pmain + phisto + plabels + paxes + ptitle + plegend + scale_fill_manual(values = getPalette(colorCount)) + pline + pfacet

## mt25
mt25.sex.rship.race = ddply(.data = idata.mt25, .variables = c("PHI_GENDER", "DEMO_DNR_RELATE", "DEMO_RACE1"), .fun = summarise, count = length(PHI_GENDER))
total2 = dcast(mt25.sex.rship.race, PHI_GENDER ~ DEMO_DNR_RELATE , sum, value.var="count")
total2 = melt(total2, id=c("PHI_GENDER"), value.name = c("count"))
total2 = plyr::rename(total2, c("variable"="DEMO_DNR_RELATE"))
total$DEMO_RACE1 = "Total"
mt25.age.eth.t = rbind(mt25.sex.rship.race, total2)
plabels = labs(x="Gender",y="Count (mt25)")
i = pmain + phisto + plabels + paxes + ptitle + plegend + scale_fill_manual(values = getPalette(colorCount)) + pline + pfacet

grid.arrange(h,i, ncol=2)
g <- arrangeGrob(h, i, ncol=2) #generates g
# ggsave("gender_race_rship_split_after_line_lt25_mt25.pdf", g, device = "pdf")

#############################################################
## agg plots
b3.dd=ddply(.data=a.d, .variables="STUDYID", .fun=transform, sum.n = length(STUDYID))
b3.d=ddply(.data = b3.dd, .variables = c("AGE_TRANSPLANT","PHI_GENDER","STUDYID"), .fun=summarise, count = length(AGE_TRANSPLANT), proportion = count/sum.n[1])
a.d = idata

## median ages
# median(a.d[(a.d$STUDYID=="ImmPort_DR19") & (a.d$PHI_GENDER=="Male"),]$AGE_TRANSPLANT)
# median(a.d[(a.d$STUDYID=="ImmPort_DR19") & (a.d$PHI_GENDER=="Female"),]$AGE_TRANSPLANT)
# 
# median(a.d[(a.d$STUDYID=="UNOS_OPTN") & (a.d$PHI_GENDER=="Male"),]$AGE_TRANSPLANT)
# median(a.d[(a.d$STUDYID=="UNOS_OPTN") & (a.d$PHI_GENDER=="Female"),]$AGE_TRANSPLANT)

## PLOT split by dataset; plot by age and PHI_GENDER using ggplot2; bars
# pmain = ggplot(b3.d, aes(AGE_TRANSPLANT, count, fill = PHI_GENDER))
b3.d = b3.d[!is.na(b3.d$PHI_GENDER),]
pmain = ggplot(b3.d, aes(AGE_TRANSPLANT, count, fill = PHI_GENDER)) ## remove unknown PHI_GENDER
phisto = geom_bar(stat = "identity", position=position_dodge())
pfacet = facet_grid(STUDYID ~., scale = "free_y") 
plabels = labs(x="Age",y="Count")
pticks = scale_x_continuous(breaks=seq(min(b3.d$AGE_TRANSPLANT),max(b3.d$AGE_TRANSPLANT),by=1))
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
pmain + phisto + pfacet + plabels + pticks + paxes + ptitle + plegend + scale_fill_brewer(palette = "Set1")
# ggsave("gender_dataset_split_after.pdf", device = "pdf")



## PLOT split by dataset; plot by age and sex using ggplot2; line plot
pmain = ggplot(b3.d, aes(AGE_TRANSPLANT, count, color = PHI_GENDER))
phisto = geom_line(size=1.5)
pfacet = facet_grid(STUDYID ~., scale = "free_y") 
plabels = labs(x="Age",y="Count")
pticks = scale_x_continuous(breaks=seq(min(b3.d$AGE_TRANSPLANT),max(b3.d$AGE_TRANSPLANT),by=1))
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
pmain + phisto + pfacet + plabels + pticks + paxes + ptitle + plegend + scale_color_brewer(palette = "Set1")
# ggsave("gender_dataset_split_after_line.pdf", device = "pdf")


