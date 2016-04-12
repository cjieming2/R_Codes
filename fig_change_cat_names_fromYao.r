filter<- read.table("C:/Users/Jieming/Desktop/Cate.Filter")
ran <- read.table("C:/Users/Jieming/Desktop/supptable1.txt",row.names=1)

name =rownames(ran)

a<- grep("bound.",name,value=FALSE)
tmp = name[a]
tmp = gsub ("^","TFM.",tmp)
tmp = gsub (".bound.extend.bed","",tmp)
name[a]=tmp


a<- grep("Proximal.bed",name,value=FALSE)
tmp = name[a]
tmp = gsub ("^","Pro.",tmp)
tmp = gsub (".Proximal.bed","",tmp)
name[a]=tmp


a<- grep("Distal.bed",name,value=FALSE)
tmp = name[a]
tmp = gsub ("^","Dis.",tmp)
tmp = gsub (".Distal.bed","",tmp)
name[a]=tmp

rownames(ran) = name


uniq.result = ran[filter[,1],]
