setwd("/Users/jiemingchen/Documents/transplantation/a_recipient/combined")

## my own library
source("/Users/jiemingchen/R_codes/jmRlib.R")

library(RImmPort)
library(DBI)
library(sqldf)
library(plyr)
library(RMySQL)
library(ggplot2)
library(gridExtra)
library(gapminder)
library(statmod)

## plot and save function
plot_age_sex = function(data,flag,range) {
  names(data)[names(data)=="Patient_Sex"]  <- "Sex"
  names(data)[names(data)=="Donation_Age"]  <- "AGE"
  names(data)[names(data)=="Patient_Race2"]  <- "RACE"
  
  ## this is a histogram, age 40 is considered at bin 39; uses un-ddply a - nope
  # pmain = ggplot(data=data, aes(x=round(data$AGE), fill=Sex))
  # phisto = geom_histogram(breaks=seq(min(round(data$AGE)),round(max(data$AGE)),by=1), position = "dodge")
  
  ## this groups each age group as a category and then takes 15 min to plot!! uses ddply b -- nope
  # pmain = ggplot(b,aes(x=SEX,y=count,fill=Organ))
  # phisto = geom_bar(stat = "identity",color="white") + facet_wrap(~AGE,nrow=1)
  
  if(flag == 1)
  {
    ############
    # this one b1 
    pmain = ggplot(data=data, aes(x=AGE, y=count, fill=Sex))
    phisto = geom_bar(stat="identity", position = position_dodge())
    plabels = labs(x="Age",y="Count")
    pticks = scale_x_continuous(breaks=seq(min(round(data$AGE)),round(max(data$AGE)),by=1))
    paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
                  axis.title.y = element_text(face = "bold",colour = "black", size = 20),
                  axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
    ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
    plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
    
    myplot = pmain + phisto + plabels + pticks + paxes + ptitle + plegend
    return(myplot)
  } else if(flag == 2)
  {
    ########
    ## this plots each x axis as an interaction of age and sex and then adjust label for age
    ## order of Sex and AGE matter - uses b
    ## plots organ with this data
    pmain = ggplot(data, aes(x = as.numeric(interaction(Sex,AGE)), y=count, fill=ORGAN))
    phisto = geom_bar(stat = "identity") 
    plabels = labs(x="Age",y="Count")
    # pticks = scale_x_continuous(breaks=range, labels=seq(min(round(data$AGE)),round(max(data$AGE)),by=1))
    #pticks = scale_x_continuous(breaks=range,labels=seq(min(round(data$AGE)),round(max(data$AGE)),by=1))
    # pticks = scale_x_continuous() ##debug
    paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
                  axis.title.y = element_text(face = "bold",colour = "black", size = 20),
                  axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
    ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
    plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
    
    pcolor = scale_fill_manual(values=c("#00BA38","#F8766D","#619CFF")) ## only for donors
    myplot = pmain + phisto + plabels + paxes + ptitle + plegend + pcolor
    return(myplot)
  } else
  {
    ########
    ## this plots each x axis as an interaction of age and sex and then adjust label for age
    ## order of Sex and AGE matter - uses b
    ## plots study with this data
    pmain = ggplot(data, aes(x = as.numeric(interaction(Sex,AGE)), y=count, fill=RACE))
    phisto = geom_bar(stat = "identity") 
    plabels = labs(x="Age",y="Count")
    # pticks = scale_x_continuous(breaks=range, labels=seq(min(round(data$AGE)),round(max(data$AGE)),by=1))
    #pticks = scale_x_continuous(breaks=range,labels=seq(min(round(data$AGE)),round(max(data$AGE)),by=1))
    # pticks = scale_x_continuous() ##debug
    paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
                  axis.title.y = element_text(face = "bold",colour = "black", size = 20),
                  axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
    ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
    plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
    
    # pcolor = scale_fill_manual(values=c("#00BA38","#F8766D","#619CFF")) ## only for donors
    myplot = pmain + phisto + plabels + paxes + ptitle + plegend
    return(myplot)
  }
}

##-----------------------------------------------------------
## input data GENDER
idata = read.table("combined_ehr_immport_unos-all-recipients-race-ethnicity-combined-gender.txt", header = T, sep = "\t")
a.d = idata
a.d$Donation_Age = round(a.d$Donation_Age)
b3.dd=ddply(.data=a.d, .variables="Data_Source", .fun=transform, sum.n = length(Data_Source))
b3.d=ddply(.data = b3.dd, .variables = c("Donation_Age","Sex","Data_Source"), .fun=summarise, count = length(Donation_Age), proportion = count/sum.n[1])

## median ages
median(a.d[(a.d$Data_Source=="ImmPort_DR19") & (a.d$Sex=="Male"),]$Donation_Age)
median(a.d[(a.d$Data_Source=="ImmPort_DR19") & (a.d$Sex=="Female"),]$Donation_Age)

median(a.d[(a.d$Data_Source=="UNOS_OPTN") & (a.d$Sex=="Male"),]$Donation_Age)
median(a.d[(a.d$Data_Source=="UNOS_OPTN") & (a.d$Sex=="Female"),]$Donation_Age)


## PLOT split by dataset; plot by age and sex using ggplot2; bars
# pmain = ggplot(b3.d, aes(Donation_Age, count, fill = Sex))
pmain = ggplot(b3.d[b3.d$Sex != "Not_Specified", ], aes(Donation_Age, count, fill = Sex)) ## remove unknown sex and age
phisto = geom_bar(stat = "identity", position=position_dodge())
pfacet = facet_grid(Data_Source ~., scale = "free_y") 
plabels = labs(x="Age",y="Count")
pticks = scale_x_continuous(breaks=seq(min(b3.d$Donation_Age),max(b3.d$Donation_Age),by=1))
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
pmain + phisto + pfacet + plabels + pticks + paxes + ptitle + plegend + scale_fill_brewer(palette = "Set1")
# ggsave("plot_combined_immport_emr_unos-jun2016_age_sex_kidney_donors_bars.pdf", device = "pdf")
# savePlot("plot_combined_immport_emr_unos-jun2016_age_sex_kidney_donors_bars-yr-1963-2007.png", type = "png")

## PLOT split by dataset; bar plots of aggregate genders
agg = ddply(.data = a.d, .variables = c("Sex","Data_Source"), .fun=summarise, count = length(Sex))
total.s = ddply(.data = a.d, .variables = c("Data_Source"), .fun=summarise, count = length(Data_Source))
agg = merge(agg, total.s, by="Data_Source"); names(agg) = c("Data_Source","Sex","Count","Total")
agg$proportion = agg$Count/agg$Total
pmain = ggplot(agg[agg$Sex != "Not_Specified", ], aes(Sex, y=proportion, fill = Sex)) ## remove unknown sex
phisto = geom_bar(stat = "identity", position=position_dodge())
pfacet = facet_grid(Data_Source ~., scale = "free_y") 
plabels = labs(x="Sex",y="Proportion")
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
pmain + phisto + pfacet + plabels + paxes + ptitle + plegend + scale_fill_brewer(palette = "Set1")
# ggsave("reci-agg_combined_immport_emr_unos-age_sex_kidney_bars.pdf", device = "pdf")


# PLOT split by dataset; plot by age and sex using ggplot2; line plot
pmain = ggplot(b3.d, aes(Donation_Age, count, color = Sex))
phisto = geom_line(size=1.5)
pfacet = facet_grid(Data_Source ~., scale = "free_y") 
plabels = labs(x="Age",y="Count")
pticks = scale_x_continuous(breaks=seq(min(b3.d$Donation_Age),max(b3.d$Donation_Age),by=1))
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
pmain + phisto + pfacet + plabels + pticks + paxes + ptitle + plegend + scale_color_brewer(palette = "Set1")
# ggsave("plot_combined_immport_emr_unos-jun2016_age_sex_kidney_donors_lines.pdf", device = "pdf")
# savePlot("plot_combined_immport_emr_unos-jun2016_age_sex_kidney_donors_lines-yr-1963-2007.png", type = "png")

## PLOT split by dataset; plot by age and sex using ggplot2; line plot proportion
pmain = ggplot(b3.d, aes(Donation_Age, proportion, color = Sex))
phisto = geom_line(size=1.5)
pfacet = facet_grid(Data_Source ~., scale = "free_y") 
plabels = labs(x="Age",y="Proportion")
pticks = scale_x_continuous(breaks=seq(min(b3.d$Donation_Age),max(b3.d$Donation_Age),by=5))
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
pmain + phisto + pfacet + plabels + pticks + paxes + ptitle + plegend + scale_color_brewer(palette = "Set1")
# ggsave("reci-unos-ehr-immport-proportion-gender-forsymp.pdf", device = "pdf")
# savePlot("plot_combined_immport_emr_unos-jun2016_age_sex_kidney_donors_lines-yr-1963-2007.png", type = "png")


## PLOT on same plot UNOS-ImmPort
sum.unos = sum(b3.d[b3.d$Data_Source=="UNOS_OPTN","count"])
sum.immport = sum(b3.d[b3.d$Data_Source=="ImmPort_DR19","count"])
pmain = ggplot(b3.d[b3.d$Data_Source=="UNOS_OPTN",], aes(Donation_Age, count/sum.unos, color = Sex))
phisto = geom_line(size=1.5, linetype=2)
phisto2 = geom_line(data=b3.d[b3.d$Data_Source=="ImmPort_DR19",], aes(Donation_Age, count/sum.immport, color = Sex))
plabels = labs(x="Age",y="Proportion")
pticks = scale_x_continuous(breaks=seq(min(b3.d$Donation_Age),max(b3.d$Donation_Age),by=5))
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
pmain + phisto + phisto2 + plabels + pticks + paxes + ptitle + plegend + scale_color_brewer(palette = "Set1")
# ggsave("plot_combined_immport_emr_unos-jun2016_age_sex_kidney_donors_lines-yr-1963-2007-unos_dotted-proportion.pdf", device = "pdf")
# ggsave("plot_combined_immport_emr_unos-jun2016_age_sex_kidney_donors_lines-yr-1987-2010-unos_dotted-proportion.pdf", device = "pdf")

# just the ImmPort DR19
pmain + phisto2 + plabels + pticks + paxes + ptitle + plegend + scale_color_brewer(palette = "Set1")
# ggsave("plot_combined_immport_age_sex_kidney_donors_lines-yr-1987-2010.pdf", device = "pdf")

## PLOT on same plot EMR-ImmPort
sum.emr = sum(b3.d[b3.d$Data_Source=="UCSF_EMR","count"])
sum.immport = sum(b3.d[b3.d$Data_Source=="ImmPort_DR19","count"])
pmain = ggplot(b3.d[b3.d$Data_Source=="UCSF_EMR",], aes(Donation_Age, count/sum.emr, color = Sex))
phisto = geom_line(size=1.5, linetype=2)
phisto2 = geom_line(data=b3.d[b3.d$Data_Source=="ImmPort_DR19",], aes(Donation_Age, count/sum.immport, color = Sex))
plabels = labs(x="Age",y="Count")
pticks = scale_x_continuous(breaks=seq(min(b3.d$Donation_Age),max(b3.d$Donation_Age),by=5))
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))
plegend = theme(legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14))
pmain + phisto + phisto2 + plabels + pticks + paxes + ptitle + plegend + scale_color_brewer(palette = "Set1")
# ggsave("plot_combined_immport_emr_unos-jun2016_age_sex_kidney_donors_lines-yr-all-emr_dotted-proportion.pdf", device = "pdf")


##---------------------------------------------
immport.b3 = b3.d[b3.d$Data_Source == "ImmPort_DR19",]
unosopt.b3 = b3.d[b3.d$Data_Source == "UNOS_OPTN",]
ucsfemr.b3 = b3.d[b3.d$Data_Source == "UCSF_EMR",]

##### fisher's test for disparity in aggregate gender
fisher.test(matrix(c(5427,4161,76502,51905),2,2), alternative="two.sided") ## ImmPort-unos
fisher.test(matrix(c(5427,4161,336,190),2,2), alternative="two.sided") ## ImmPort-EHR

prop.test(matrix(c(5427,76502,4161,51905),2,2)) ## ImmPort-unos
prop.test(matrix(c(5427,336,4161,190),2,2)) ## ImmPort-EHR



#### merge immport and unos #### 
## MALE
a = merge(unosopt.b3[unosopt.b3$Sex=="Male",], immport.b3[immport.b3$Sex=="Male",], by="Donation_Age", all.x=TRUE, all.y=TRUE)
a = a[,c("Donation_Age","count.x","count.y")]
a[is.na(a)] = 0

## chi sq test, expected = unos
z = a[a$count.x>20 & a$count.y>20,]
chisq.test(z[,2:3]) ## x can be count data, but y needs to be proportion NOT counts

## KS test
prob.x = a$count.x/sum(a$count.x)
prob.y = a$count.y/sum(a$count.y)
ks.test(prob.x,prob.y, exact=TRUE, alternative="two.sided")

## qqplot
qqplot(prob.x, prob.y)
abline(0,1, col="red", lty=3)

#
## FEMALE
a = merge(unosopt.b3[unosopt.b3$Sex=="Female",], immport.b3[immport.b3$Sex=="Female",], by="Donation_Age", all.x=TRUE, all.y=TRUE)
a = a[,c("Donation_Age","count.x","count.y")]
a[is.na(a)] = 0

## chi sq test, expected = unos
z = a[a$count.x>20 & a$count.y>20,]
prob.y = z$count.y/sum(z$count.y)
chisq.test(z[,2:3]) ## x can be count data, but y needs to be proportion NOT counts

## KS test
prob.x = a$count.x/sum(a$count.x)
prob.y = a$count.y/sum(a$count.y)
ks.test(prob.x,prob.y, exact=TRUE, alternative="two.sided")

## qqplot
qqplot(prob.x, prob.y)
abline(0,1, col="red", lty=3)

#
## fishers between male and female
a = merge(immport.b3[immport.b3$Sex=="Male",], immport.b3[immport.b3$Sex=="Female",], by="Donation_Age", all.x=TRUE, all.y=TRUE)
a = a[,c("Donation_Age","count.x","count.y")]
names(a) = c("Donation_Age","male.immport","female.immport")
a = merge(a, unosopt.b3[unosopt.b3$Sex=="Male",], by="Donation_Age", all.x=TRUE, all.y=TRUE)
a = merge(a, unosopt.b3[unosopt.b3$Sex=="Female",], by="Donation_Age", all.x=TRUE, all.y=TRUE)
a = a[,c("Donation_Age","male.immport","female.immport","count.x","count.y")]
names(a) = c("Donation_Age","male.immport","female.immport","male.unos","female.unos")
a[is.na(a)] = 0

fish <- function(x) {
  aa=
  bb=cbind(aa$p.value,aa$estimate)
  return(bb)
}

ft.r.p = as.data.frame(apply(a,1,function(x) fisher.test(matrix(x[2:5],2,2), alternative="two.sided")$p.value))
ft.r.e = as.data.frame(apply(a,1,function(x) fisher.test(matrix(x[2:5],2,2), alternative="two.sided")$estimate))
ft.r = cbind(a$Donation_Age,ft.r.p,ft.r.e); names(ft.r) = c("Donation_Age","fish.p","fish.OR_est")
ft.r$fish.p.bon = as.data.frame(apply(as.matrix(ft.r$fish.p),1,function(x) min(x*nrow(ft.r),1)) )

#### merge emr and immport #### 
## MALE
a = merge(ucsfemr.b3[ucsfemr.b3$Sex=="Male",], immport.b3[immport.b3$Sex=="Male",], by="Donation_Age", all.x=TRUE, all.y=TRUE)
a = a[,c("Donation_Age","count.x","count.y")]
a[is.na(a)] = 0

## chi sq test, expected = unos
z = a[a$count.x>0 & a$count.y>0,]
prob.y = z$count.y/sum(z$count.y)
chisq.test(z[,2:3]) ## x can be count data, but y needs to be proportion NOT counts

## KS test
prob.x = a$count.x/sum(a$count.x)
prob.y = a$count.y/sum(a$count.y)
ks.test(prob.x,prob.y, exact=TRUE, alternative="two.sided")

## qqplot
qqplot(prob.x, prob.y)
abline(0,1, col="red", lty=3)

#
## FEMALE
a = merge(ucsfemr.b3[ucsfemr.b3$Sex=="Female",], immport.b3[immport.b3$Sex=="Female",], by="Donation_Age", all.x=TRUE, all.y=TRUE)
a = a[,c("Donation_Age","count.x","count.y")]
a[is.na(a)] = 0

## chi sq test, expected = unos
z = a
prob.y = z$count.y/sum(z$count.y)
chisq.test(z$count.x, prob.y) ## x can be count data, but y needs to be proportion NOT counts

## KS test
prob.x = a$count.x/sum(a$count.x)
prob.y = a$count.y/sum(a$count.y)
ks.test(prob.x,prob.y, exact=TRUE, alternative="two.sided")

## qqplot
qqplot(prob.x, prob.y)
abline(0,1, col="red", lty=3)

#
## fishers between male and female
a = merge(ucsfemr.b3[ucsfemr.b3$Sex=="Male",], immport.b3[immport.b3$Sex=="Female",], by="Donation_Age", all.x=TRUE, all.y=TRUE)
a = a[,c("Donation_Age","count.x","count.y")]
names(a) = c("Donation_Age","male.immport","female.immport")
a = merge(a, ucsfemr.b3[ucsfemr.b3$Sex=="Male",], by="Donation_Age", all.x=TRUE, all.y=TRUE)
a = merge(a, ucsfemr.b3[ucsfemr.b3$Sex=="Female",], by="Donation_Age", all.x=TRUE, all.y=TRUE)
a = a[,c("Donation_Age","male.immport","female.immport","count.x","count.y")]
names(a) = c("Donation_Age","male.immport","female.immport","male.emr","female.emr")
a[is.na(a)] = 0

fish <- function(x) {
  aa=
    bb=cbind(aa$p.value,aa$estimate)
  return(bb)
}

ft.r.p = as.data.frame(apply(a,1,function(x) fisher.test(matrix(x[2:5],2,2), alternative="two.sided")$p.value))
ft.r.e = as.data.frame(apply(a,1,function(x) fisher.test(matrix(x[2:5],2,2), alternative="two.sided")$estimate))
ft.r = cbind(a$Donation_Age,ft.r.p,ft.r.e); names(ft.r) = c("Donation_Age","fish.p","fish.OR_est")
ft.r$fish.p.bon = as.data.frame(apply(as.matrix(ft.r$fish.p),1,function(x) min(x*nrow(ft.r),1)) )  ## this file combined_emr_immport_unos-year-1963-2007.txt has no EMR data




#### merge emr and unos #### 
## MALE
a = merge(ucsfemr.b3[ucsfemr.b3$Sex=="Male",], unosopt.b3[unosopt.b3$Sex=="Male",], by="Donation_Age", all.x=TRUE, all.y=TRUE)
a = a[,c("Donation_Age","count.x","count.y")]
a[is.na(a)] = 0

## chi sq test, expected = unos
# z = a[a$count.x>20 & a$count.y>20,]
z=a
prob.y = z$count.y/sum(z$count.y)
chisq.test(z$count.x, prob.y) ## x can be count data, but y needs to be proportion NOT counts

## KS test
prob.x = a$count.x/sum(a$count.x)
prob.y = a$count.y/sum(a$count.y)
ks.test(prob.x,prob.y, exact=TRUE, alternative="two.sided")

## qqplot
qqplot(prob.x, prob.y)
abline(0,1, col="red", lty=3)

#
## FEMALE
a = merge(ucsfemr.b3[ucsfemr.b3$Sex=="Female",], unosopt.b3[unosopt.b3$Sex=="Female",], by="Donation_Age", all.x=TRUE, all.y=TRUE)
a = a[,c("Donation_Age","count.x","count.y")]
a[is.na(a)] = 0

## chi sq test, expected = unos
z = a
prob.y = z$count.y/sum(z$count.y)
chisq.test(z$count.x, prob.y) ## x can be count data, but y needs to be proportion NOT counts

## KS test
prob.x = a$count.x/sum(a$count.x)
prob.y = a$count.y/sum(a$count.y)
ks.test(prob.x,prob.y, exact=TRUE, alternative="two.sided")

## qqplot
qqplot(prob.x, prob.y)
abline(0,1, col="red", lty=3)

#
## fishers between male and female
a = merge(ucsfemr.b3[ucsfemr.b3$Sex=="Male",], ucsfemr.b3[ucsfemr.b3$Sex=="Female",], by="Donation_Age", all.x=TRUE, all.y=TRUE)
a = a[,c("Donation_Age","count.x","count.y")]
names(a) = c("Donation_Age","male.emr","female.emr")
a = merge(a, unosopt.b3[unosopt.b3$Sex=="Male",], by="Donation_Age", all.x=TRUE, all.y=TRUE)
a = merge(a, unosopt.b3[unosopt.b3$Sex=="Female",], by="Donation_Age", all.x=TRUE, all.y=TRUE)
a = a[,c("Donation_Age","male.emr","female.emr","count.x","count.y")]
names(a) = c("Donation_Age","male.emr","female.emr","male.unos","female.unos")
a[is.na(a)] = 0

fish <- function(x) {
  aa=
    bb=cbind(aa$p.value,aa$estimate)
  return(bb)
}

ft.r.p = as.data.frame(apply(a,1,function(x) fisher.test(matrix(x[2:5],2,2), alternative="two.sided")$p.value))
ft.r.e = as.data.frame(apply(a,1,function(x) fisher.test(matrix(x[2:5],2,2), alternative="two.sided")$estimate))
ft.r = cbind(a$Donation_Age,ft.r.p,ft.r.e); names(ft.r) = c("Donation_Age","fish.p","fish.OR_est")
ft.r$fish.p.bon = as.data.frame(apply(as.matrix(ft.r$fish.p),1,function(x) min(x*nrow(ft.r),1)) ) ## this file combined_emr_immport_unos-year-1963-2007.txt has no EMR data
