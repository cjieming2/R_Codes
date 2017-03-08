setwd('/Users/jiemingchen/Documents/varimed')
library(RCurl)
library(rjson)

# Read sample UUIDs.
uuids <- read.table(file="uuid2sample2cancertype.txt", sep="\t")

# Convert to character vector.
uuids <- as.vector(t(uuids))

# Query TCGA's UUID to barcode Web Service.
resp <- getURL("https://tcga-data.nci.nih.gov/uuid/uuidws/mapping/json/uuid/batch", customrequest="POST", httpheader=c("Content-Type: text/plain"), postfields=paste(uuids, collapse=","))

# Extract mappings from response.
mappings <- fromJSON(resp)$uuidMapping

# Extract patient barcode from sample barcode.
mappings <- lapply(mappings, function (mapping) {
  mapping$barcode <- substr(mapping$barcode, 0, 12)
  return(mapping)
})
