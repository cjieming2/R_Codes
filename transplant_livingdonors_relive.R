setwd("/Users/jiemingchen/Documents/transplantation/a_donor/immport/studyfiles/final")

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
## input data PRE and POST TRANSPLANT
pre.t = read.delim("merged-relive01_02_03_immport2tsv-noControls_pretransplant-exam_kidney.txt", header = T, sep = "\t")
post.t = read.delim("merged-relive01_02_03_immport2tsv-noControls_posttransplant-exam_kidney.txt", header = T, sep = "\t")
demo = read.delim("merged-relive01_02_03_immport2tsv-noControls_demo-social_kidney.txt", header = T, sep = "\t")

## Height; pre.t. EVAL_HGT = in cm
# compare pre and post tx height - they should be almost the same!
pre.t.height = select(pre.t, Sub_Org_Accession, EVAL_HGT)
pre.t.height = pre.t.height[complete.cases(pre.t.height),]
post.t.height = select(post.t, Sub_Org_Accession, HEIGHT_COMBINED)
post.t.height = post.t.height[complete.cases(post.t.height),]
pre.post.height = inner_join(pre.t.height, post.t.height, by="Sub_Org_Accession")
pmain = ggplot(pre.post.height, aes(EVAL_HGT, HEIGHT_COMBINED, label = Sub_Org_Accession))
phisto = geom_point()
plabels = labs(x="pre.tx.height",y="post.tx.height")
ptext = geom_text(aes(label=ifelse(abs(HEIGHT_COMBINED-EVAL_HGT)>20,as.character(Sub_Org_Accession),'')), vjust = 0, hjust = 0, size = 5)
pmain + phisto + plabels + paxes + ptitle + plegend + ptext
# ggsave("pre.post.tx.height.pdf", device = "pdf")


## Weight
## egrep "Sub_|SUB124608|SUB127063|SUB127094|SUB127321" merged-relive01_02_03_immport2tsv-noControls_posttransplant-exam_kidney.txt | ftranspose | grep -n '' | egrep "WEIGHT|WGT|Sub_|GFR_DO_|ALLEXAMDONESAMEDAY|CHS_TODAY" | sed 's/:/\t/g' | cut -f2 | egrep "GFR_DO|ALLEXAMDONESAMEDAY|CHS_TODAY"
post.t.bmi.weight.dates = select(post.t, Sub_Org_Accession, GFR_DO_GFR1_MIN, GFR_DO_GFR2_MIN, GFR_DO_GFR3_MIN, GFR_DO_GFR4_MIN, GFR_DO_GFR1002_MIN, GFR_DO_GFR2002_MIN, GFR_DO_GFR3002_MIN, GFR_DO_GFR4002_MIN, GFR_DO_GFR1003_MIN, GFR_DO_GFR2003_MIN, GFR_DO_GFR3003_MIN, GFR_DO_GFR4003_MIN, GFR_DO_GFR1004_MIN, GFR_DO_GFR2004_MIN, GFR_DO_GFR3004_MIN, GFR_DO_GFR4004_MIN, GFR_DO_GFR1005_MIN, GFR_DO_GFR2005_MIN, GFR_DO_GFR3005_MIN, GFR_DO_GFR4005_MIN, GFR_DO_GFR1006_MIN, GFR_DO_GFR2006_MIN, GFR_DO_GFR3006_MIN, GFR_DO_GFR4006_MIN, GFR_DO_GFR1007_MIN, GFR_DO_GFR2007_MIN, GFR_DO_GFR3007_MIN, GFR_DO_GFR4007_MIN, GFR_DO_GFR1008_MIN, GFR_DO_GFR2008_MIN, GFR_DO_GFR3008_MIN, GFR_DO_GFR4008_MIN, WEIGHTDATEDONE_MIN, CHS_TODAY_MIN)
post.t.bmi.weight.dates = cbind.data.frame(post.t.bmi.weight.dates,1,2)
names(post.t.bmi.weight.dates) = c("Sub_Org_Accession", "GFR_WGT1", "GFR_WGT2", "GFR_WGT3", "GFR_WGT4", "GFR_WGT1002", "GFR_WGT2002", "GFR_WGT3002", "GFR_WGT4002", "GFR_WGT1003", "GFR_WGT2003", "GFR_WGT3003", "GFR_WGT4003", "GFR_WGT1004", "GFR_WGT2004", "GFR_WGT3004", "GFR_WGT4004", "GFR_WGT1005", "GFR_WGT2005", "GFR_WGT3005", "GFR_WGT4005", "GFR_WGT1006", "GFR_WGT2006", "GFR_WGT3006", "GFR_WGT4006", "GFR_WGT1007", "GFR_WGT2007", "GFR_WGT3007", "GFR_WGT4007", "GFR_WGT1008", "GFR_WGT2008", "GFR_WGT3008", "GFR_WGT4008", "WEIGHTVALUE", "CHS_WEIGHT", "WEIGHT_COMBINED", "GFRWEIGHT")
post.t.bmi.weight = select(post.t, Sub_Org_Accession, GFR_WGT1, GFR_WGT2, GFR_WGT3, GFR_WGT4, GFR_WGT1002, GFR_WGT2002, GFR_WGT3002, GFR_WGT4002, GFR_WGT1003, GFR_WGT2003, GFR_WGT3003, GFR_WGT4003, GFR_WGT1004, GFR_WGT2004, GFR_WGT3004, GFR_WGT4004, GFR_WGT1005, GFR_WGT2005, GFR_WGT3005, GFR_WGT4005, GFR_WGT1006, GFR_WGT2006, GFR_WGT3006, GFR_WGT4006, GFR_WGT1007, GFR_WGT2007, GFR_WGT3007, GFR_WGT4007, GFR_WGT1008, GFR_WGT2008, GFR_WGT3008, GFR_WGT4008, WEIGHTVALUE, CHS_WEIGHT, WEIGHT_LB, WEIGHT_KG, GFRWEIGHT)
post.t.bmi.weight = mutate(post.t.bmi.weight, WEIGHT_COMBINED = ifelse(is.na(WEIGHT_KG), ifelse(is.na(WEIGHT_LB), NA, WEIGHT_LB*0.45), WEIGHT_KG) )
post.t.bmi.weight = select(post.t.bmi.weight, Sub_Org_Accession, GFR_WGT1, GFR_WGT2, GFR_WGT3, GFR_WGT4, GFR_WGT1002, GFR_WGT2002, GFR_WGT3002, GFR_WGT4002, GFR_WGT1003, GFR_WGT2003, GFR_WGT3003, GFR_WGT4003, GFR_WGT1004, GFR_WGT2004, GFR_WGT3004, GFR_WGT4004, GFR_WGT1005, GFR_WGT2005, GFR_WGT3005, GFR_WGT4005, GFR_WGT1006, GFR_WGT2006, GFR_WGT3006, GFR_WGT4006, GFR_WGT1007, GFR_WGT2007, GFR_WGT3007, GFR_WGT4007, GFR_WGT1008, GFR_WGT2008, GFR_WGT3008, GFR_WGT4008, WEIGHTVALUE, CHS_WEIGHT, WEIGHT_COMBINED, GFRWEIGHT)
post.t.bmi.weight[post.t.bmi.weight == 0] = NA

post.t.bmi.weight.dates.1 = melt(post.t.bmi.weight.dates, id.vars = "Sub_Org_Accession", variable.name= "datetype", value.name = "date")
post.t.bmi.weight.1 = melt(post.t.bmi.weight, id.vars = "Sub_Org_Accession", variable.name= "weighttype", value.name = "weight"); 

post.t.weight = inner_join(post.t.bmi.weight.dates.1, post.t.bmi.weight.1, by = c("Sub_Org_Accession" = "Sub_Org_Accession", "datetype" = "weighttype"))
post.t.weight$weight[post.t.weight$weight > 700] = NA
post.t.weight$weight[post.t.weight$weight < 3] = NA # look out for min weight for babies: 2kg
post.t.weight = post.t.weight[!is.na(post.t.weight$weight) & post.t.weight$weight != 0,]

# ptest = post.t.weight[post.t.weight$Sub_Org_Accession == "SUB123554" | post.t.weight$Sub_Org_Accession == "SUB123555" | post.t.weight$Sub_Org_Accession == "SUB123556" | post.t.weight$Sub_Org_Accession == "SUB123557" | post.t.weight$Sub_Org_Accession == "SUB123558" | post.t.weight$Sub_Org_Accession == "SUB123559" | post.t.weight$Sub_Org_Accession == "SUB123560", ]
# ptest = ptest[!is.na(ptest$weight) & ptest$weight != 0,]
# pmain = ggplot(ptest, aes(datetype, weight))
# phisto = geom_line(aes(colour=Sub_Org_Accession, group=Sub_Org_Accession))
# phisto2 = geom_point(aes(colour=Sub_Org_Accession), size = 3)
## track trajectory of post tx weight change
## realize there are missing data
pmain = ggplot(post.t.weight, aes(datetype, weight))
phisto = geom_line(aes(colour=Sub_Org_Accession, group=Sub_Org_Accession))
phisto2 = geom_point(aes(colour=Sub_Org_Accession), size = 3)
plabels = labs(x="datetype",y="weight (kg)")
## set1 usually has only 8 colors
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colorCount = length(unique(ptest$Sub_Org_Accession))
## "manually" add colors instead of scale_fill_brewer fixed brwer colors
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
pmain + phisto + phisto2 + plabels + paxes + 
  scale_fill_manual(values = getPalette(colorCount)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(legend.position="none")
# ggsave("post.tx.weight.trajectory.pdf", device = "pdf")


## missing data for which weight
nc = ncol(post.t.bmi.weight) - 1
nonmissingrate = as.data.frame((nc - apply(is.na(post.t.bmi.weight),1,sum))); names(nonmissingrate) = "nonmissingrate"
nmcount = ddply(.data = nonmissingrate, .variables = c("nonmissingrate"), .fun = summarise, count = length(nonmissingrate))
pmain = ggplot(nmcount, aes(nonmissingrate,count))
phisto = geom_bar(stat = "identity")
# phisto = geom_histogram(breaks=seq(min(nonmissingrate$nonmissingrate), max(nonmissingrate)+1, by =1), col="red", fill="green")
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
pticks = scale_x_continuous(breaks=seq(min(nonmissingrate), max(nonmissingrate)+1, by =1))
pmain + phisto + paxes + pticks + geom_text(data=nmcount,aes(nonmissingrate,y=count,label=count),vjust=-1, size = 5)
# pmain + phisto + paxes + pticks + coord_cartesian(xlim=c(1, max(nonmissingrate)+1), ylim=c(0,500))
# ggsave("post.tx.weight.missingdata.pdf", device = "pdf")



## BMI https://www.nhlbi.nih.gov/health/educational/lose_wt/BMI/bmicalc.htm
# Underweight = <18.5
# Normal weight = 18.5–24.9
# Overweight = 25–29.9
# Obesity = BMI of 30 or greater

pre.tt = inner_join(pre.t, demo, by="Sub_Org_Accession")
pre.bmi = select(pre.tt, Sub_Org_Accession, EVAL_WGT, EVAL_HGT, AGE_TRANSPLANT)
pre.bmi = mutate(pre.bmi, bmi = EVAL_WGT / (EVAL_HGT/100 * EVAL_HGT/100) )
pre.bmi = pre.bmi[complete.cases(pre.bmi),]
pmain = ggplot(pre.bmi, aes(bmi, fill=AGE_TRANSPLANT))
# phisto = geom_histogram(breaks=seq(18,65,by=1), aes(fill=AGE_TRANSPLANT))
phisto = geom_bar(stat="bin")
pline1 = geom_vline(xintercept=18.5, color="red")
pline2 = geom_vline(xintercept=25, color="red")
pline3 = geom_vline(xintercept=30, color="red")
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))

pmain + phisto + pline1 + pline2 + pline3 + paxes + ptitle 
# ggsave("pretransplant_bmi_histogram.pdf", device = "pdf")



## 2D plot btwn bmi and age
pmain = ggplot(pre.bmi, aes(AGE_TRANSPLANT,bmi))
phisto = geom_point()
pline1 = geom_vline(xintercept=20, color="red")
pline2 = geom_hline(yintercept=25, color="red")
pline3 = geom_hline(yintercept=30, color="red")
pline4 = geom_hline(yintercept=18.5, color="red")
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))

pmain + phisto + pline1 + pline2 + pline3 + pline4 + paxes + ptitle 
# ggsave("pretransplant_bmi_age.pdf", device = "pdf")







### recheck post.t weights and heights


##-----------------------------------------------------------##-----------------------------------------------------------##-----------------------------------------------------------##-----------------------------------------------------------##-----------------------------------------------------------##-----------------------------------------------------------##-----------------------------------------------------------##-----------------------------------------------------------##-----------------------------------------------------------##-----------------------------------------------------------##-----------------------------------------------------------##---------------------------------------------------------
##-----------------------------------------------------------##-----------------------------------------------------------##-----------------------------------------------------------
## input data DEMO SOCIAL
idata = read.delim("merged-relive01_02_03_immport2tsv-noControls_demo-social_kidney.txt", header = T, sep = "\t")
idata=mutate(idata, AGE_TRANSPLANT=round(AGE_TRANSPLANT)) 

##
## add in recipient information
recip.data = read.delim("merged-relive01_02_03_immport2tsv-noControls_recipient_kidney.txt", header = T, sep = "\t")

## extract recipient relationship
recip.rship = select(recip.data, Sub_Org_Accession, DEMO_DNR_RELATE, DEMO_OTH_NONBIO, DEMO_OTH_BIO)
idata.reciprship = inner_join(idata, recip.rship, by = "Sub_Org_Accession")

## split at age 25
idata.lt25 = idata.reciprship[idata.reciprship$AGE_TRANSPLANT <= 25,]
idata.mt25 = idata.reciprship[idata.reciprship$AGE_TRANSPLANT > 25,]

##-----------------------------------------------------------
## relationships
age.rship = ddply(.data = idata.reciprship, .variables = c("AGE_TRANSPLANT", "DEMO_DNR_RELATE"), .fun=summarise, count = length(AGE_TRANSPLANT))
pmain = ggplot(age.rship, aes(AGE_TRANSPLANT, count, fill = DEMO_DNR_RELATE))
phisto = geom_bar(stat = "identity")
plabels = labs(x="Age",y="Count")
pticks = scale_x_continuous(breaks=seq(min(age.rship$AGE_TRANSPLANT),max(age.rship$AGE_TRANSPLANT),by=1))
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
pline = geom_vline(xintercept=25, color="red")
## set1 usually has only 8 colors
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
## how many colors do we need
colorCount = length(unique(age.rship$DEMO_DNR_RELATE))
## "manually" add colors instead of scale_fill_brewer fixed brwer colors
pmain + phisto + plabels + pticks + paxes + ptitle + plegend + scale_fill_manual(values = getPalette(colorCount)) + pline
# ggsave("age_reciprship_bar.pdf", device = "pdf")

#
# line plot
pmain = ggplot(age.rship, aes(AGE_TRANSPLANT, count, color = DEMO_DNR_RELATE))
phisto = geom_line(size=1.5)
pmain + phisto + plabels + pticks + paxes + ptitle + plegend + scale_color_manual(values = getPalette(colorCount)) + pline
# ggsave("age_reciprship_line.pdf", device = "pdf")

#
# split by gender
age.rship.gender = ddply(.data = idata.reciprship, .variables = c("AGE_TRANSPLANT", "PHI_GENDER", "DEMO_DNR_RELATE"), .fun=summarise, count = length(AGE_TRANSPLANT))
pmain = ggplot(age.rship.gender, aes(AGE_TRANSPLANT, count, fill = DEMO_DNR_RELATE))
phisto = geom_bar(stat = "identity")
pfacet = facet_grid(PHI_GENDER ~., scale = "free_y")
pmain + phisto + plabels + pticks + paxes + ptitle + plegend + scale_fill_manual(values = getPalette(colorCount)) + pline + pfacet
# ggsave("age_gender_split_reciprship_bar.pdf", device = "pdf")

pmain = ggplot(age.rship.gender, aes(AGE_TRANSPLANT, count, color = DEMO_DNR_RELATE))
phisto = geom_line(size=1.5)
pfacet = facet_grid(PHI_GENDER ~., scale = "free_y")
pmain + phisto + plabels + pticks + paxes + ptitle + plegend + scale_color_manual(values = getPalette(colorCount)) + pline + pfacet
# ggsave("age_gender_split_reciprship_line.pdf", device = "pdf")


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


