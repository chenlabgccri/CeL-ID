## Script for cell line identification. Requires DP and/or Freq files for 934 as well as the query cell line.
## Before running this script please make sure that both the files should have chromosome represented in similar style.
## chr1, chr2, ... chry or 1, 2, ...y

library(data.table)
library(corrplot)
## Read DP and FREQ files for 934 cell lines
dp <- fread("CCLE_DP.txt", sep  = "\t")
#freq <- fread("CCLE_FREQ.txt", sep  = "\t")

## Read DP and/or FREQ file for the query cell line
qDp <- read.delim("Query_DP.txt", sep  = "\t")
#qFreq<- read.delim("Query_FREQ.txt", sep  = "\t")

## Select all the variants present in the 934 cell lines and merging
d.merge <- merge(dp, qDp, by = "Position", all.x = T)
#d.merge <- merge(freq, qFreq, by = "Position", all.x = T)

d.merge[is.na(d.merge)] = 0

## Merging based on variants in the query cell line
#d.merge <- merge(dp, newdp, by = "Position", all.y = T)
#d.merge[is.na(d.merge)] = 0

## correlation calculation between query and all 934 cell lines
s =data.frame(cor(d.merge$Sample1, d.merge[,-1]))

## Best match
s = s[,-935]
corrplot(as.matrix(sort(s, decreasing = TRUE)[1:5]))
colnames(as.matrix(sort(s, decreasing = TRUE)[1]))
