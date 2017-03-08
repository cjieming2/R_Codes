setwd("/Users/jiemingchen/Documents/transplantation/a_donor")

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
## input data
idata = read.table("combined_emr_immport_unos-all.txt", header = T, sep = "\t")
a.d = idata
b3.d=ddply(.data = a.d, .variables = c("Donation_Age","Sex","Data_Source"), .fun=summarise, count = length(Donation_Age))

## PLOT split by dataset; plot by age and sex using ggplot2; bars
pmain = ggplot(b3.d, aes(Donation_Age, count, fill = Sex))
# pmain = ggplot(b3.d[b3.d$Sex != "Not_Specified", ], aes(Donation_Age, count, fill = Sex)) ## remove unknown sex
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
# savePlot("plot_combined_immport_emr_unos-jun2016_age_sex_kidney_donors_bars.png", type = "png")

## PLOT split by dataset; plot by age and sex using ggplot2; line plot
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
# savePlot("plot_combined_immport_emr_unos-jun2016_age_sex_kidney_donors_lines.png", type = "png")


##---------------------------------------------
immport.b3 = b3.d[b3.d$Data_Source == "ImmPort_DR19",]
unosopt.b3 = b3.d[b3.d$Data_Source == "UNOS_OPTN",]
ucsfemr.b3 = b3.d[b3.d$Data_Source == "UCSF_EMR",]

#### merge immport and unos #### 
## MALE
a = merge(unosopt.b3[unosopt.b3$Sex=="Male",], immport.b3[immport.b3$Sex=="Male",], by="Donation_Age", all.x=TRUE, all.y=TRUE)
a = a[,c("Donation_Age","count.x","count.y")]
a[is.na(a)] = 0

## chi sq test, expected = unos
z = a[a$count.x>20 & a$count.y>20,]
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
a = merge(unosopt.b3[unosopt.b3$Sex=="Female",], immport.b3[immport.b3$Sex=="Female",], by="Donation_Age", all.x=TRUE, all.y=TRUE)
a = a[,c("Donation_Age","count.x","count.y")]
a[is.na(a)] = 0

## chi sq test, expected = unos
z = a[a$count.x>20 & a$count.y>20,]
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
ft.r$fish.p.bon = as.data.frame(apply(as.matrix(ft.r$fish.p),1,function(x) min(x*nrow(ft.r),1)) )

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
ft.r$fish.p.bon = as.data.frame(apply(as.matrix(ft.r$fish.p),1,function(x) min(x*nrow(ft.r),1)) )
