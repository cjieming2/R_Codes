#
#    Please read PubMed E-utility Usage Guidelines here: http://www.ncbi.nlm.nih.gov/books/NBK25497/
#
#    Check for updates and read more about this script at: http://rpsychologist.com  
#
#    Packages used: "XML", "Rcurl" and "plyr" 


# -----------Run these lines first---------------------
### You have to set working directory to were scripts are located ###
setwd("/path/do/directory/")

# Script source
source("PubMedTrend.R")

# -----------------------------------------------------

# Example Search Query
query <- c("cbt" = "cognitive behav* psychotherap*[tiab] OR cognitive behav* therap*[tiab] OR cognitive therap*[tiab] OR cognitive psychotherap*[tiab] OR behav* therap*[tiab] OR behav* psychotherapy[tiab]",
           "pdt" = "psychodynamic therap*[tiab] OR psychodynamic psychotherap*[tiab] OR psychoanalytic therap*[tiab] OR psychoanalytic psychotherap*[tiab]")

# Run the script
df <- PubMedTrend(query)

#  Run this to get the total hits for each query in 'query'
#     Specify if you want to get the whole 'query' or just the 'names'
#     enter no argument  if you want booth 'query' and 'name'. 
df.hits <- PubTotalHits()
df.hits <- PubTotalHits("query")
df.hits <- PubTotalHits("names")


#######################
#### EXAMPLE PLOTS ####
#######################

# Note:
# ------
# These plots wont work if you only have 1 query. 

library("ggplot2")
library("directlabels")

### AREA PLOT ###
ggplot(df, aes(year, relative, group=.id, fill=.id)) + 
  geom_area() +
  opts(title=paste("Area Plot of PubMed Publications per Year\nfor", paste(names(query), collapse = ", "))) +
  xlab("year") +
  ylab("Publications per 1 million PubMed articles") +
  scale_fill_brewer()

### LINE PLOTS ###

# RAW
ggplot(df, aes(year, relative, group=.id, color=.id)) + 
  geom_line(show_guide=F) + 
  xlab("Publication year") +
  ylab("Publications per 1 million PubMed articles") +
  opts(title = paste("Pubmed hits for", paste(names(query), collapse = ", ")))

# SMOOTHED
p <- ggplot(df, aes(year, relative, group=.id, color=.id)) + 
  geom_line(alpha = I(7/10), color="grey", show_guide=F) +
  stat_smooth(size=2, span=0.3, se=F, show_guide=F) + 
  xlab("Publication year") +
  ylab("Publications per 1 million PubMed articles") +
  opts(title = paste("Pubmed hits (smoothed) for", paste(names(query), collapse = ", "))) +
  xlim(1950,2020)
direct.label(p, "last.bumpup")