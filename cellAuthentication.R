#    Demo code for Cel-ID functions.
#     
##   Yidong Chen, Tabrez Mohammad Oct 2019, GCCRI/UTHSCSA

if (!requireNamespace("data.table", quietly = TRUE)){
  install.packages("data.table")}
library(data.table)

if (!requireNamespace("corrplot", quietly = TRUE)){
  install.packages("corrplot")}
library(corrplot)

source( "Cel-ID_functions.R")

## Read DP and FREQ files for 934 cell lines
ccle_dp <- read.delim("CCLE_DP.txt", sep  = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
ccle_freq <- read.delim("CCLE_FREQ.txt", sep  = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

## Read DP and/or FREQ file for the query cell line
query_dp <- read.delim("Query_DP.txt", sep  = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
query_freq<- read.delim("Query_FREQ.txt", sep  = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

# find the best match, sdx codes the sorted match order based on correlation coefficient
rt = CCLE_Identification( ccle_freq, ccle_dp, query_freq, query_dp )

# return with correlation coefficient, sorted order index, and top/bottom 5 correlation matrix.
rho = rt[[1]]
sdx = rt[[2]]
corMatrixTop5Bottom5 = rt[[3]]

# plot correlation matrix.
corrplot( corMatrixTop5Bottom5, method = "pie", type = "lower")


##
# find possible contaminator
# create a fake mixtures (70% targeted cell, 30% contaminator)
matched_ccleID = sdx[[1]]
mixtureID = 10
numPos = nrow( ccle_freq )
p1 = rnorm( numPos, 0.70, 0.05 );
p2 = rnorm( numPos, 0.30, 0.02 );

# make sure freq is within [0,100], and DP is greater than and equal to 0 after adding noise
mixFreq = data.frame( ccle_freq[, matched_ccleID] * p1 + ccle_freq[, mixtureID] * p2 )
rownames(mixFreq) <- rownames(ccle_freq)
colnames(mixFreq) <- c("mix")
mixFreq[mixFreq[,1] < 0, 1] = 0
mixFreq[mixFreq[,1] > 100, 1] = 100

mixDP = data.frame( trunc( ccle_dp[, matched_ccleID] * p1  +  ccle_dp[, mixtureID] * p2 ) )
rownames(mixDP) <- rownames(ccle_freq)
colnames(mixDP) <- c("mix")
mixFreq[mixDP[,1] < 0, 1] = 0

# test. 
rt = CCLE_MixtureEstimate( ccle_freq, ccle_dp, mixFreq, mixDP, sdx[1] )


