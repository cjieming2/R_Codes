## ========================================================================== ##
## T2D-GENES Project 1 phenotype_GWAS cleaning
## First check on phenotype data
## Check for variables for which data is available for each cohort.
## For those that are available, obtain the summary statistics for each, including the number of samples with missing information
## ========================================================================== ##

# snowwhite:/net/snowwhite/home/xlsim/t2d-genes/Project1/phenotypes_GWAS> ls
# Ashkenazi  Jackson_Heart_Study  Korea  LOLIPOP  METSIM  San_Antonio  Singapore_East_Asian  Singapore_South_Asian  Starr_County  Wake_Forest

root.read.dir <- "~/t2d-genes/Project1/phenotypes_GWAS/src"
root.write.dir <- "~/t2d-genes/Project1/phenotypes_GWAS/phenotype"
popname.dir <- c("Ashkenazi","Jackson_Heart_Study","Korea","LOLIPOP","METSIM","San_Antonio","Singapore_East_Asian","Singapore_South_Asian","Starr_County","Wake_Forest")
pheno.file <- c("NA","T2DGENES_project1_JHS_phenotypes_selected_sequenced_samples.csv","East_Asian_Korea_Phenotype_20111007.csv","pheno.txt","metsim_pheno.csv","NA","phenotype/Singapore-EastAsian-SDCSSP2-phenotype-variables-dcc-01082011.txt","phenotype/Singapore-SouthAsian-SINDI-phenotype-variables-dcc-01082011.txt","HS_pheno.txt","WFU_GWAS_T2Dgenes.phe")


date <- "25oct2011"
varlist <- c("FAMID","STUDYID","FATHER","MOTHER","SEX","AGE","T2D","AOD","FAST_GLU","FAST_INS","FAST_CPEP","2H_GLU","2H_INS","2H_CPEP","HBA1C","GAD","CREATININE","ADIPONECTIN","LEPTIN","CHOL","LDL","HDL","TG","HEIGHT","WEIGHT","BMI","HIPC","WAISTC","SBP","DBP","DIABMEDS","BPMEDS","LIPIDMEDS","HORMONES")


## Function to get everyone onto the same set of variables and headers from the colnames provided by cohorts #
fill.table.function <- function(file.in)
	{
	dat.new <- data.frame(array(data=NA,dim=c(nrow(dat),length(varlist)),dimnames=list(NULL,varlist)),stringsAsFactors=F)
	missing.col <- NULL
	for (j in 1:ncol(dat.new))
		{
		id <- which(colnames(file.in)%in%colnames(dat.new)[j])
		if (length(id)!=0) dat.new[,j] <- file.in[,id]
		if (length(id)==0) missing.col <- c(missing.col,colnames(dat.new)[j])
		}
	return(list(dat.new, missing.col))
	}
# Function to check if there are any duplicated FAMID, FATHER or MOTHER #
family.check.function <- function(file.in)
  {
  fam.dup.id <- which(!is.na(file.in[,"FAMID"]) & duplicated(file.in[,"FAMID"])==T)
  fam.dup <- which(file.in[,"FAMID"]%in%file.in[fam.dup.id,"FAMID"])
  father.dup.id <- which(!is.na(file.in[,"FATHER"]) & duplicated(file.in[,"FATHER"])==T)
  father.dup <- which(file.in[,"FATHER"]%in%file.in[father.dup.id,"FATHER"])
  mother.dup.id <- which(!is.na(file.in[,"MOTHER"]) & duplicated(file.in[,"MOTHER"])==T)
  mother.dup <- which(file.in[,"MOTHER"]%in%file.in[mother.dup.id,"MOTHER"])
  fam.id <- unique(c(fam.dup,father.dup,mother.dup))
  if (length(fam.id)!=0)
    {
    family <- file.in[fam.id,]
    family.sort <- family[order(family[,"FAMID"]),]
    write.table(family.sort[,1:4], file=paste(root.write.dir,"/",popname.dir[i],"_related_individuals_spotted_",date,".txt",sep=""), quote=F, row.names=F, sep="\t")
    }
  }
# Function to check for variables that were provided but had NA values for all samples #
var_avail_function <- function(file.in)
  {
  tmp1 <- NULL
  for (j in 1:ncol(file.in))
    {
    tmp <- ifelse(all(is.na(file.in[,j])==T),0,1)
    tmp1 <- c(tmp1, tmp)
    }
  var.n.avail <- names(file.in)[which(tmp1==0)]
  return(var.n.avail)
  }
## Function to compute the summary stats for each of the cohorts, all and by case control status ##
#  continuous traits: N, mean, median (sd) [min, max] NA
#  binary traits: N (%) males, yes.
summarize.function <- function(file.in, pop)
  {
  cts <- c("AGE","AOD","FAST_GLU","FAST_INS","FAST_CPEP","2H_GLU","2H_INS","2H_CPEP","HBA1C","GAD","CREATININE","ADIPONECTIN","LEPTIN","CHOL","LDL","HDL","TG","HEIGHT","WEIGHT","BMI","HIPC","WAISTC","SBP","DBP")
  bin <- c("SEX","T2D","DIABMEDS","BPMEDS","LIPIDMEDS","HORMONES")
  row.id <- which(rownames(summary.info)%in%pop)
  j=1
  while (j <= ncol(file.in))
    {
    vari <- colnames(file.in)[j]
    if (vari%in%cts)
      {
      if (all(is.na(file.in[,j]))!=T)
	      {
	      col.id <- which(colnames(summary.info)%in%vari)
	      ntotal <- length(file.in[,j])
	      n <- length(which(!is.na(file.in[,j])))
	      mean1 <- mean(file.in[,j], na.rm=T)
	      median1 <- median(file.in[,j], na.rm=T)
	      sd1 <- median(file.in[,j], na.rm=T)
	      min1 <- min(file.in[,j], na.rm=T)
	      max1 <- max(file.in[,j], na.rm=T)
	      na <- ntotal - n
	      out <- paste("N=",n,",mean=",formatC(mean1,format="f",digits=2),",median=",formatC(median1,format="f",digits=2),",sd=",formatC(sd1,format="f",digits=2),",min=",formatC(min1,format="f",digits=2),",max=",formatC(max1,format="f",digits=2),",NA=",na,sep="")
	      summary.info[row.id,col.id] <- out
	      #print(out)
	      #print(c(row.id,col.id))
	      }
      }
    if (vari%in%bin)
      {
      if (all(is.na(file.in[,j]))!=T)
	      {
	      col.id <- which(colnames(summary.info)%in%vari)
	      ntotal <- length(file.in[,j])
	      n <- length(which(!is.na(file.in[,j])))
	      n1 <- length(which(file.in[,j]==1))
	      n1.percent <- n1/n
	      na <- ntotal - n
	      out <- paste("N=",ntotal,",1=",n1,",1%=",formatC(n1.percent,format="f",digits=2),",0=",n-n1,",0%=",formatC(1-n1.percent,format="f",digits=2),",NA=",na,sep="")
	      summary.info[row.id,col.id] <- out
	      #print(out)
	      #print(c(row.id,col.id))
	      }
      }
    j=j+1
    #print(j)
    }
  return(summary.info)
  }
  

## Tables of information to fill in ##
cohort.info <- data.frame(study=popname.dir, var_not_provided=rep(NA,times=10), var_not_avail=rep(NA,times=10),stringsAsFactors=F)
summary.info <- as.data.frame(matrix(NA,nrow=30,ncol=30,dimnames=list(sort(c(popname.dir,paste(popname.dir,"cases",sep="."),paste(popname.dir,"controls",sep="."))),varlist[5:length(varlist)])),stringsAsFactors=F)


# Jackson Heart Study #
# ------------------- #
i=2
filename.tmp <- paste(root.read.dir,popname.dir[i],pheno.file[i],sep="/")
dat <- read.table(filename.tmp,header=T,sep=",",stringsAsFactors=F)
dat[1:3,]
# Check, format header and print missing variables #
dat.new <- fill.table.function(file.in=dat)[[1]]
missing.col <- fill.table.function(file.in=dat)[[2]]
missing.col
names(dat.new) <- varlist
if (length(missing.col)!=0) cohort.info[i,2] <- paste(missing.col,collapse=", ")
# Check for variables with NA for all samples #
var.avail <- var_avail_function(dat.new)
if (length(var.avail)!=0) cohort.info[i,3] <- paste(var.avail,collapse=", ")
# Check for families #
family.check <- family.check.function(dat.new)
# Check T2D case/control ascertainment and if it agrees with the T2D variable
# Fill in the summary characteristics #
summary.info <- summarize.function(dat.new, popname.dir[i])
dat.new.cases <- dat.new[which(dat.new[,"T2D"]==1),]
summary.info <- summarize.function(dat.new.cases, paste(popname.dir[i],"cases",sep="."))
dat.new.controls <- dat.new[which(dat.new[,"T2D"]==0),]
summary.info <- summarize.function(dat.new.controls, paste(popname.dir[i],"controls",sep="."))

# Korea #
# ----- #
i=3
filename.tmp <- paste(root.read.dir,popname.dir[i],pheno.file[i],sep="/")
dat <- read.table(filename.tmp,header=T,sep=",",stringsAsFactors=F)
dat[1:3,]
# Check, format header and print missing variables #
dat.new <- fill.table.function(file.in=dat)[[1]]
missing.col <- fill.table.function(file.in=dat)[[2]]
missing.col
names(dat)[which(names(dat)%in%"FAST_ING")] <- "FAST_INS"
names(dat)[which(names(dat)%in%"LIPMEDS")] <- "LIPIDMEDS"
dat.new <- fill.table.function(file.in=dat)[[1]]
missing.col <- fill.table.function(file.in=dat)[[2]]
if (length(missing.col)!=0) cohort.info[i,2] <- paste(missing.col,collapse=", ")
# Check for variables with NA for all samples #
var.avail <- var_avail_function(dat.new)
if (length(var.avail)!=0) cohort.info[i,3] <- paste(var.avail,collapse=", ")
names(dat.new) <- varlist
# Check for families #
family.check <- family.check.function(dat.new)
# Check T2D case/control ascertainment and if it agrees with the T2D variable
dat.new$t2d <- NA
t2d1.id <- which(dat.new[,"DIABMEDS"]==1)
t2d2.id <- which(dat.new[,"FAST_GLU"]>=7)
t2d3.id <- which(dat.new[,"2H_GLU"]>=11.1)
t2d4.id <- which(dat.new[,"AOD"]>=40)
dat.new[unique(c(t2d1.id,t2d2.id,t2d3.id,t2d4.id)),"t2d"] <- 1
not2d.id <- which((dat.new[,"DIABMEDS"]==0 | is.na(dat.new[,"DIABMEDS"])) & (dat.new[,"FAST_GLU"]<5.6 & dat.new[,"2H_GLU"]<7.8))
dat.new[not2d.id,"t2d"] <- 0
dat.new[which(dat.new[,"T2D"]!=dat.new[,"t2d"]),]
dat.new <- dat.new[,1:34]
# Fill in the summary characteristics #
summary.info <- summarize.function(dat.new, popname.dir[i])
dat.new.cases <- dat.new[which(dat.new[,"T2D"]==1),]
summary.info <- summarize.function(dat.new.cases, paste(popname.dir[i],"cases",sep="."))
dat.new.controls <- dat.new[which(dat.new[,"T2D"]==0),]
summary.info <- summarize.function(dat.new.controls, paste(popname.dir[i],"controls",sep="."))

# LOLIPOP #
# ------- #
i=4
filename.tmp <- paste(root.read.dir,popname.dir[i],pheno.file[i],sep="/")
dat <- read.table(filename.tmp,header=T,sep="\t",stringsAsFactors=F) # special case of tab instead of csv
dat[1:3,]
# Check, format header and print missing variables #
names(dat) <- toupper(names(dat))
dat.new <- fill.table.function(file.in=dat)[[1]]
missing.col <- fill.table.function(file.in=dat)[[2]]
missing.col
paste(missing.col,collapse=", ")
names(dat)[which(names(dat)%in%"SAMPLE")] <- "STUDYID"
names(dat)[which(names(dat)%in%"DM")] <- "T2D"
names(dat)[which(names(dat)%in%"GLUCOSE")] <- "FAST_GLU"
names(dat)[which(names(dat)%in%"INSULIN")] <- "FAST_INS"
names(dat)[which(names(dat)%in%"TRIG")] <- "TG"
names(dat)[which(names(dat)%in%"HIP")] <- "HIPC"
names(dat)[which(names(dat)%in%"WAIST")] <- "WAISTC"
names(dat)[which(names(dat)%in%"DMHX")] <- "DIABMEDS"
names(dat)[which(names(dat)%in%"HTHX")] <- "BPMEDS"
dat.new <- fill.table.function(file.in=dat)[[1]]
missing.col <- fill.table.function(file.in=dat)[[2]]
if (length(missing.col)!=0) cohort.info[i,2] <- paste(missing.col,collapse=", ")
# Check for variables with NA for all samples #
var.avail <- var_avail_function(dat.new)
if (length(var.avail)!=0) cohort.info[i,3] <- paste(var.avail,collapse=", ")
names(dat.new) <- varlist
# Check for families #
family.check <- family.check.function(dat.new)
# Fill in the summary characteristics #
summary.info <- summarize.function(dat.new, popname.dir[i])
dat.new.cases <- dat.new[which(dat.new[,"T2D"]==1),]
summary.info <- summarize.function(dat.new.cases, paste(popname.dir[i],"cases",sep="."))
dat.new.controls <- dat.new[which(dat.new[,"T2D"]==0),]
summary.info <- summarize.function(dat.new.controls, paste(popname.dir[i],"controls",sep="."))

# METSIM #
# ------ #
i=5
filename.tmp <- paste(root.read.dir,popname.dir[i],pheno.file[i],sep="/")
dat <- read.table(filename.tmp,header=T,sep=",",stringsAsFactors=F) # special case of tab instead of csv
dat[1:3,]
# Check, format header and print missing variables #
dat.new <- fill.table.function(file.in=dat)[[1]]
missing.col <- fill.table.function(file.in=dat)[[2]]
missing.col
if (length(missing.col)!=0) cohort.info[i,2] <- paste(missing.col,collapse=", ")
# Check for variables with NA for all samples #
var.avail <- var_avail_function(dat.new)
if (length(var.avail)!=0) cohort.info[i,3] <- paste(var.avail,collapse=", ")
names(dat.new) <- varlist
# Check for families #
family.check <- family.check.function(dat.new)
# Fill in the summary characteristics #
summary.info <- summarize.function(dat.new, popname.dir[i])
dat.new.cases <- dat.new[which(dat.new[,"T2D"]==1),]
summary.info <- summarize.function(dat.new.cases, paste(popname.dir[i],"cases",sep="."))
dat.new.controls <- dat.new[which(dat.new[,"T2D"]==0),]
summary.info <- summarize.function(dat.new.controls, paste(popname.dir[i],"controls",sep="."))

# Singapore_East_Asian #
# -------------------- #
i=7
filename.tmp <- paste(root.read.dir,popname.dir[i],pheno.file[i],sep="/")
dat <- read.table(filename.tmp,header=T,sep=",",stringsAsFactors=F) # special case of tab instead of csv
dat[1:3,]
# Check, format header and print missing variables #
dat.new <- fill.table.function(file.in=dat)[[1]]
missing.col <- fill.table.function(file.in=dat)[[2]]
missing.col
names(dat)[which(names(dat)%in%"HIP")] <- "HIPC"
dat.new <- fill.table.function(file.in=dat)[[1]]
missing.col <- fill.table.function(file.in=dat)[[2]]
missing.col
if (length(missing.col)!=0) cohort.info[i,2] <- paste(missing.col,collapse=", ")
# Check for variables with NA for all samples #
var.avail <- var_avail_function(dat.new)
if (length(var.avail)!=0) cohort.info[i,3] <- paste(var.avail,collapse=", ")
names(dat.new) <- varlist
# Check for families #
family.check <- family.check.function(dat.new)
# Fill in the summary characteristics #
summary.info <- summarize.function(dat.new, popname.dir[i])
dat.new.cases <- dat.new[which(dat.new[,"T2D"]==1),]
summary.info <- summarize.function(dat.new.cases, paste(popname.dir[i],"cases",sep="."))
dat.new.controls <- dat.new[which(dat.new[,"T2D"]==0),]
summary.info <- summarize.function(dat.new.controls, paste(popname.dir[i],"controls",sep="."))

# Singapore_South_Asian #
# --------------------- #
i=8
filename.tmp <- paste(root.read.dir,popname.dir[i],pheno.file[i],sep="/")
dat <- read.table(filename.tmp,header=T,sep=",",stringsAsFactors=F) # special case of tab instead of csv
dat[1:3,]
# Check, format header and print missing variables #
dat.new <- fill.table.function(file.in=dat)[[1]]
missing.col <- fill.table.function(file.in=dat)[[2]]
missing.col
names(dat)[which(names(dat)%in%"HIP")] <- "HIPC"
names(dat)[which(names(dat)%in%"H_GLU")] <- "X2H_GLU"
names(dat)[which(names(dat)%in%"H_INS")] <- "X2H_INS"
names(dat)[which(names(dat)%in%"H_CPEP")] <- "X2H_CPEP"
dat.new <- fill.table.function(file.in=dat)[[1]]
missing.col <- fill.table.function(file.in=dat)[[2]]
missing.col
if (length(missing.col)!=0) cohort.info[i,2] <- paste(missing.col,collapse=", ")
# Check for variables with NA for all samples #
var.avail <- var_avail_function(dat.new)
if (length(var.avail)!=0) cohort.info[i,3] <- paste(var.avail,collapse=", ")
names(dat.new) <- varlist
# Check for families #
family.check <- family.check.function(dat.new)
# Fill in the summary characteristics #
summary.info <- summarize.function(dat.new, popname.dir[i])
dat.new.cases <- dat.new[which(dat.new[,"T2D"]==1),]
summary.info <- summarize.function(dat.new.cases, paste(popname.dir[i],"cases",sep="."))
dat.new.controls <- dat.new[which(dat.new[,"T2D"]==0),]
summary.info <- summarize.function(dat.new.controls, paste(popname.dir[i],"controls",sep="."))

# Starr Country #
# ------------- #
i=9
filename.tmp <- paste(root.read.dir,popname.dir[i],pheno.file[i],sep="/")
dat <- read.table(filename.tmp,header=T,sep=",",stringsAsFactors=F) # special case of tab instead of csv
dat[1:3,]
# Check, format header and print missing variables #
dat.new <- fill.table.function(file.in=dat)[[1]]
missing.col <- fill.table.function(file.in=dat)[[2]]
missing.col
if (length(missing.col)!=0) cohort.info[i,2] <- paste(missing.col,collapse=", ")
# Check for variables with NA for all samples #
var.avail <- var_avail_function(dat.new)
if (length(var.avail)!=0) cohort.info[i,3] <- paste(var.avail,collapse=", ")
names(dat.new) <- varlist
# Check for families #
family.check <- family.check.function(dat.new)
# Fill in the summary characteristics #
summary.info <- summarize.function(dat.new, popname.dir[i])
dat.new.cases <- dat.new[which(dat.new[,"T2D"]==1),]
summary.info <- summarize.function(dat.new.cases, paste(popname.dir[i],"cases",sep="."))
dat.new.controls <- dat.new[which(dat.new[,"T2D"]==0),]
summary.info <- summarize.function(dat.new.controls, paste(popname.dir[i],"controls",sep="."))


# Wake Forest #
# ----------- #
i=10
filename.tmp <- paste(root.read.dir,popname.dir[i],pheno.file[i],sep="/")
dat <- read.table(filename.tmp,header=T,sep=",",stringsAsFactors=F) # special case of tab instead of csv
dat[1:3,]
# Check, format header and print missing variables #
dat.new <- fill.table.function(file.in=dat)[[1]]
missing.col <- fill.table.function(file.in=dat)[[2]]
missing.col
if (length(missing.col)!=0) cohort.info[i,2] <- paste(missing.col,collapse=", ")
# Check for variables with NA for all samples #
var.avail <- var_avail_function(dat.new)
if (length(var.avail)!=0) cohort.info[i,3] <- paste(var.avail,collapse=", ")
names(dat.new) <- varlist
# Check for families #
family.check <- family.check.function(dat.new)
# Fill in the summary characteristics #
summary.info <- summarize.function(dat.new, popname.dir[i])
dat.new.cases <- dat.new[which(dat.new[,"T2D"]==1),]
summary.info <- summarize.function(dat.new.cases, paste(popname.dir[i],"cases",sep="."))
dat.new.controls <- dat.new[which(dat.new[,"T2D"]==0),]
summary.info <- summarize.function(dat.new.controls, paste(popname.dir[i],"controls",sep="."))


write.table(cohort.info, file=paste(root.write.dir,"/missing_variable_by_cohort_",date,".txt",sep=""), quote=F, row.names=F, col.names=T, sep="\t")
summary.info.new <- data.frame(cohort=rownames(summary.info),summary.info,stringsAsFactors=F)
write.table(summary.info.new, file=paste(root.write.dir,"/summary_statistics_by_cohort_and_T2D_",date,".txt",sep=""), quote=F, row.names=F, col.names=T, sep="\t")

