setwd("C:/Shared/scripts-R_perl_shell_macros/R codes/PubMedTrend");

# the scripts are downloaded from 
# http://rpsychologist.com/an-r-script-to-automatically-look-at-pubmed-citation-counts-by-year-of-publication/

# script source
source("PubMedTrend.R");

# query sample
# query <- c("cbt"= "cognitive behav* psychotherap*[tiab] OR cognitive behav* therap*[tiab]",
#            "pdt" = "psychodynamic therap*[tiab] OR psychodynamic psychotherap*[tiab]",
#            "psychoanalytic" = "psychoanalytic therap*[tiab] OR psychoanalytic psychoterap*[tiab]",
#            "ssri" = "selective serotonin reuptake inhibitor*[tiab]",
#            "mindfulness" = "mindfulness[tiab]")

# query <- c("tpr" = "tpr [tw]","ank" = "ankyrin [tw]","armadillo" = "armadillo[tw]");
## note that you need to change the years manually in the function script [PPDAT] vs [DP]
# query <- c("tpr" = "tpr [tw] OR tetratricopeptide [tw]","ank" = "ank [tw] OR ankyrin [tw]",
#            "armadillo" = "armadillo [tw]", "lrr"="lrr [tw] OR leucine rich repeat [tw]",
#            "sh3" = "sh3 [tw]");

# query <- c("KEGG" = "KEGG [tw]","NetPath" = "NetPath [tw]","Reactome" = "Reactome [tw]","SignaLink" = "SignaLink [tw]");
query <- c("YY Teo" = "yik ying teo");

df <- PubMedTrend(query, yrMax = 2020);

# show number of hits
pubhits <- PubTotalHits();

# plot
x11()
library(ggplot2)
library(directlabels)
ggplot(df, aes(year, relative, group=.id, fill=.id)) +
  geom_area() +
  opts(title=paste("Area Plot of PubMed Publications per Year\nfor", paste(names(query), collapse = ", "))) +
  xlab("year") +
  ylab("Publications per 1 million PubMed articles") +
  scale_fill_brewer()

### LINE PLOTS ###
x11()
# RAW
ggplot(df, aes(year, relative, group=.id, color=.id)) +
  geom_line(show_guide=F) +
  xlab("Publication year") +
  ylab("Publications per 1 million PubMed articles") +
  opts(title = paste("Pubmed hits for", paste(names(query), collapse = ", ")))

x11()
# SMOOTHED
p <- ggplot(df, aes(year, relative, group=.id, color=.id)) +
  geom_line(alpha = I(7/10), color="grey", show_guide=F) +
  stat_smooth(size=2, span=0.3, se=F, show_guide=F) +
  xlab("Publication year") +
  ylab("Publications per 1 million PubMed articles") +
  opts(title = paste("Pubmed hits (smoothed) for", paste(names(query), collapse = ", "))) +
  xlim(1950,2020)
direct.label(p, "last.bumpup")