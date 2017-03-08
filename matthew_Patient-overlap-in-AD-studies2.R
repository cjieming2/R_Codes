library(RImmPort) 
library(DBI) 
library(sqldf) 
library(plyr)
library(RMySQL)
library(dplyr)


# provide MySQL database parameters
mysql_conn <- dbConnect(MySQL(), user='lab_user', password='big_data', dbname='proj_study_ALLSTUDIES', host='exxport-mysql.celiaramktyo.us-west-2.rds.amazonaws.com')

# set data source as ImmPort MySQL database
setImmPortDataSource(mysql_conn)

#########################################

# list of AD studies
AD_studies <- c('SDY2', 'SDY4', 'SDY5', 'SDY6', 'SDY7', 'SDY8', 'SDY9', 'SDY10', 'SDY13', 'SDY14')

#list of Allergic Rhinitis studies
AR_studies <- c('SDY1', 'SDY223')

#list of Food Allergy studies

FA_studies <- c('SDY218', 'SDY660')

#list of Asthma studies
# stud(LB)

Asthma_studies <- c('SDY210', 'SDY211', 'SDY545')

##########################################



##########################################
#Get Unique Subject IDs for AD studies using for loop and vector of lists called AD_uniqueIDs

AD_uniqueIDs <- c()

# Sets domain for demographics
domain_name <- "Demographics"

# Populates AD_uniqueIDs a list of values from each study
for (i in 1:length(AD_studies)) {
  demographics <- getDomainDataOfStudies(domain_name, AD_studies[i])
  if (length(demographics) > 0)
    names(demographics)
  AD_uniqueIDs[i] <- list(demographics$dm_df["USUBJID"])
}

# defines vector for length of studies
AD_studyLengths <- c()

# Gives lengths of the studys
for (i in 1:length(AD_studies)) {
  AD_studyLengths <- c(AD_studyLengths, length(AD_uniqueIDs[[i]][[1]]))
  
  # reports the length of study i in the list of vectors
  # length(AD_uniqueIDs[[i]][[1]])
  
}

# generates data frame with Study Name
AD_Names_Lengths <- data.frame(AD_studies, stringsAsFactors = FALSE)
# adds new column with lengths of study
AD_Names_Lengths$Patients <- AD_studyLengths
# Outputs a data.frame with the name of the study, number of patients per study
AD_Names_Lengths


# output specific value in AD_uniqueIDs
AD_uniqueIDs[[1]][[1]]
# AD_uniqueIDs[[x]][[1]][[y]]
# where x = study number and y = patient number

###### Generating patient overlap lists ########################

AD_overlapIDs <- c()
k <- 1
# length(AD_uniqueIDs[[2]][[1]])

# testing output of just 1 cell
# AD_uniqueIDs[[1]][[1]][[1]]

for (i in 1:length(AD_uniqueIDs)) {
  for (j in 1:length(AD_uniqueIDs)) {
    
    # adds list of overlaps between study i and j to AD_overlapIDs
    AD_overlapIDs[k] <- list(intersect(AD_uniqueIDs[[i]][[1]], AD_uniqueIDs[[j]][[1]]))
    k <- k + 1 # counts to 100
    
  } 
  
}

AD_overlapIDs

lengths_overlaps <- c()

# gives lengths of overlaps in AD_overlapIDs
for (i in 1:length(AD_overlapIDs)){
  lengths_overlaps[i] <- c(length(AD_overlapIDs[[i]]))
}
 lengths_overlaps

 
# makes a data frame for numbers of overlaps
x <- 1:10
y <- 1:10

# empty data frame 10 x 10
overlaps <- data.frame(x)
for (i in 1:10) {
  overlaps[i] <- x
}
overlaps

# populates data frame with number of overlaps from
k <- 1
for (i in 1:10) {
  for (j in 1:10) {
    overlaps[i, j] <- lengths_overlaps[k]
    k <- k + 1
  }
}
overlaps

# renames the columns and row names of overlaps to study name
names(overlaps) <- AD_studies
row.names(overlaps) <- AD_studies
overlaps # finished product over overlaps!!

#####################################################

# makes a data frame of studies with total numbers in study

# starts with size and labeling of overlaps

totals_of_studies <- overlaps
# AD_studyLengths is the lengths of the studies in order

for(i in 1:10) {
  for(j in 1:10) {
    
    totals_of_studies[i, j] <- AD_studyLengths[j]
    
  }
}

totals_of_studies

#######################################################

# generates number of unique patients

uniquePatients <- 0
for (i in 1:10) {
  uniquePatients <- uniquePatients + overlaps[i, i] * 2 - sum(overlaps[i:10, i])
}
uniquePatients # 1751 unique patients!

# generates total numbers of subjects (sum of diagonal)
total_subjects <- 0
for (i in 1:10) {
  total_subjects <- total_subjects + overlaps[i, i]
}
total_subjects # 2282 subjects collected

uniquePatients/total_subjects # 76.7% of patients unique


########################################################

# Square Pie function
squarePie <- function(pct, col="black", col.grid="#e0e0e0", col.border="black", main="") {
  
  if (pct > 100) {
    pct <- 100
    warning("Percentage value, pct, should be an integer between 0 and 100")
  } else if (pct < 0) {
    pct <- 0
    warning("Percentage value, pct, should be an integer between 0 and 100.")
  }
  
  # Round to nearest integer
  pct <- round(pct)
  
  # x- and y-coordinates of rows and columns
  x_row <- 1:10
  y_col <- 1:10
  
  # put together full coordinate vectors
  x <- rep(x_row, 10)
  y <- rep(y_col, each=10)
  
  # set colors
  fill_col <- c(rep(col, pct), rep("#ffffff", 100-pct))
  
  # plot
  plot(0, 0, type="n", xlab="", ylab="", main=main, xlim=c(0, 11), ylim=c(0, 10.5), asp=1, bty="n", axes=FALSE)
  symbols (x, y, asp=1, squares=rep(1,100), inches=FALSE, add=TRUE, bg=fill_col, fg=col.grid, lwd=0.5)
  rect(.5, .5, 10.5, 10.5, lwd=2, border=col.border)
  
}

###############################################

# Generating 100 square Pis

percentage_overlap <- overlaps / totals_of_studies * 100
percentage_overlap

par(mfrow = c(10, 10), mar = c(0, 0, 0, 0))

for(i in 1:10) {
  for(j in 1:10) {
    border_color = "#B4B9BF"
    squarePie(percentage_overlap[i, j], col.grid=NA, col="#007CBE", col.border=border_color)
  }
}

# messinga round with colors

par(mfrow = c(10, 10), mar = c(0, 0, 0, 0), x)

for(i in 1:10) {
  for(j in 1:10) {
    if (i!=j) {
      border_color = "#B4B9BF"
      squarePie(percentage_overlap[i, j], col.grid=NA, col="#007CBE", col.border=border_color)
    } else {
      border_color = "#000000"
      squarePie(percentage_overlap[i, j], col.grid=NA, col="#e0e0e0", col.border=border_color)
    }
    if (overlaps[i, j] > 0) {
      text(5.5, 5.5, overlaps[i, j], font=2, cex=1.5, family = "Roboto")
    }
  }
}


