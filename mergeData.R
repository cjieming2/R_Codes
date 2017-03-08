setwd('/Users/jiemingchen/Documents/transplantation/a_donor/immport/studyfiles/final/src')

## for some reason read.table does not detect same number of columns for some rows of the files
data1 = read.delim("relive_01_data_jcedited_decoded.txt_", header = T, sep = "\t", stringsAsFactors = FALSE, na.strings = "")
data2 = read.delim("relive_02_data_jcedited_decoded.txt_", header = T, sep = "\t", stringsAsFactors = FALSE, na.strings = "")

data12 = merge(data1, data2, all = TRUE, sort = TRUE)
write.table(data12, file="merged-relive01_02.txt", sep = "\t", quote = FALSE)


data3 = read.delim("relive_03_data_jcedited_decoded.txt_", header = T, sep = "\t", stringsAsFactors = FALSE, na.strings = "")
data123 = merge(data12, data3, all = TRUE, sort = TRUE)
write.table(data123, file="merged-relive01_02_03.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## jm is a copy of immport_DR19_DM_immport2tsv_donors_9962subjects_8studies.mergeDup.organ.txt, with USUBJID --> Sub_Org_Accession
data4 = read.delim("jm", header = T, sep = "\t", stringsAsFactors = FALSE, na.strings = "")
data1234 = merge(data123, data4, all = TRUE, sort = TRUE)
write.table(data1234, file="merged-relive01_02_03_immport2tsv.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## ----########----------

setwd('/Users/jiemingchen/Documents/transplantation/a_donor/immport/studyfiles/final/src')