#   CCLE_Identification( ccle_freq, ccle_dp, query_freq, query_dp )
#
# Args:
#     ccle_freq:         SNP frequency of all cell-lines in CCLE collection
#     ccle_dp:           SNP depth of coverage of CCLE collection
#     query_freq:        SNP frequency of query cell-line
#     query_dp:          SNP depth of coverage of query cell-line
#     
##   Yidong Chen, Tabrez Mohammad Oct 2019, GCCRI/UTHSCSA

CCLE_Identification <- function( ccle_freq, ccle_dp, query_freq, query_dp ) {

  ## Select all the variants present in the 934 cell lines and merging
  freq_merged <- merge(ccle_freq, query_freq, by = "row.names", all = FALSE )
  dp_merged   <- merge(ccle_dp, query_dp, by = "row.names", all = FALSE )
  
  rownames(freq_merged) <- freq_merged$Row.names;
  freq_merged = freq_merged[,-1]
  rownames(dp_merged) <- dp_merged$Row.names;
  dp_merged = dp_merged[,-1]
  numRows = nrow( freq_merged )
  numSamples = ncol( freq_merged )
  cat( paste( "Total of ", numRows, " variants selected" ), "\n" )

  # initialization (note: last sample are query sample)
  rho <- rep( 0, numSamples-1 )
  jcount <- rep( 0, numSamples-1 )
  snpCounts <- rep(0, numSamples-1 )
  dpThr <- 10
  
  tmp1 = freq_merged > 0;
  tmp2 = dp_merged >= dpThr;

  for( i in 1:(numSamples-1) ) {
    snpCounts[i] = sum( tmp1[,i] & tmp2[,i] )
    idx = which( (tmp1[,i] | tmp1[,numSamples] ) & (tmp2[,i] & tmp2[,numSamples]) )
    jcount[i] = length(idx)
    if( jcount[i] <= 1 ) {
      rho[i] = 0;
    } else {
      rho[i] = cor( freq_merged[idx,i], freq_merged[idx, numSamples] )
    }
  }
  
  # output each cell-line and it best and second to best cell-line name and
  # correlation coefficient/number of matching mutations.
  sdx = sort( rho, decreasing = TRUE, index.return = TRUE )$ix
  
  # Within CCLE RNA-seq samples, rho distribution is (mu, sigma, p0.001) = 0.470, 0.047, 0.61
  # see T. Mohammad, et al, BMC Genomics, 2019
  txtStr = c("1st", "2nd", "3rd")
  cell_names = colnames(freq_merged)
  colnames(freq_merged) <- c( cell_names[-numSamples], "Query sample")
  for( k in 1:3 ) {
    p = pnorm( rho[sdx[k]], mean = 0.470, sd = 0.047, lower.tail = FALSE );
    if( p > 0.001 ) {
      cat( paste(txtStr[k], " match cell-line is ", cell_names[sdx[k]], 
                            "with rho = ", sprintf( "%.3f", rho[sdx[k]]), 
                            " p = ",  sprintf( "%.3f", p) ) )
    } else {
      cat( paste(txtStr[k], " match cell-line is ", cell_names[sdx[k]], 
                            "with rho = ", sprintf( "%.3f", rho[sdx[k]]), 
                            " p = ",  sprintf( "%.3e", p) ) )
    }
    if( rho[sdx[k]] > 0.61 && k == 1 ) {
      cat( "   <<< matched cell-line (> L0.001)" )
    }
    cat( "\n")
  }
  
  m = freq_merged[,c(numSamples, sdx[1:5], sdx[(numSamples-5):(numSamples-1)])]
  for( i in 1:6 ) {
    for( j in 1:6 ) {
      if( m[i,j] == 0 ) { m[i, j] = NA}
    }
  }
  mx = cor( m, use = "pairwise.complete.obs" )
  
  # return correlation/sort order index as list.
  return( list( rho, sdx, mx ) )
}

#   CCLE_MixtureEstimate( ccle_freq, ccle_dp, query_freq, query_dp, matched_ccleID )
#
# Args:
#     ccle_freq:         SNP frequency of all cell-lines in CCLE collection
#     ccle_dp:           SNP depth of coverage of CCLE collection
#     query_freq:        SNP frequency of query cell-line
#     query_dp:          SNP depth of coverage of query cell-line
#     matched_ccleID:    The target CCLE sample ID, to be tested in combination with all other CCLEs    
#
##   Yidong Chen, Tabrez Mohammad Oct 2019, GCCRI/UTHSCSA
CCLE_MixtureEstimate <- function( ccle_freq, ccle_dp, query_freq, query_dp, matched_ccleID ) {
  
  # parameter initialization. We require the depth of coverage to be at least 10 fold
  dpThr = 10
  
  cellNames = colnames(ccle_freq)
  numSamples = ncol( ccle_freq )
  predict = rep( 0, numSamples )
  predict2 = rep( 0, numSamples )
  ts = rep( 0, numSamples )
  pval = rep( 1, numSamples )
  
  # merge all samples to the common SNPs
  freq_merged <- merge(ccle_freq, query_freq, by = "row.names", all = FALSE )
  dp_merged   <- merge(ccle_dp, query_dp, by = "row.names", all = FALSE )
  rownames(freq_merged) <- freq_merged$Row.names;
  freq_merged = freq_merged[,-1]
  rownames(dp_merged) <- dp_merged$Row.names;
  dp_merged = dp_merged[,-1]
  
  # get the numSamples updated (last one is the query sample)
  numRows = nrow( freq_merged )
  numSamples = ncol( freq_merged )
  
  # loop through all ccle samples (except the matched_ccleID, and any other cell-lines 
  # with high correlation coefficient with the matched_ccleID
  tsSaved = 0
  for (guessedMix  in  1:(numSamples - 1)) {
    if (guessedMix != matched_ccleID) {
      x1 = freq_merged[, matched_ccleID]
      x2 = freq_merged[, guessedMix]
      y = freq_merged[, numSamples]
      
      # Within CCLE RNA-seq samples, rho distribution is (mu, sigma, p0.001) = 0.470, 0.047, 0.61
      # see T. Mohammad, et al, BMC Genomics, 2019
      rho = cor(x1, x2)
      if (rho < 0.61) {
        dp1 = dp_merged[, matched_ccleID]
        dp2 = dp_merged[, guessedMix]
        idx = which(dp1 >= dpThr & dp2 >= dpThr)
        
        x1 = x1[idx]
        x2 = x2[idx]
        y = y[idx]
        
        # we select 200 most positive differential SNPs, middle 200, and 200 negative differential SNPs
        aa = x1 - x2
        adx = sort(aa, index.return = TRUE)$ix
        
        numMuts = length(aa)
        midPoint = trunc(numMuts / 2)
        idx = c(adx[1:200], adx[(midPoint - 100):(midPoint + 100)], adx[(numMuts - 200):numMuts])
        
        # build the data frame for lm()
        testData = data.frame(x1 = x1[idx], x2 = x2[idx], y = as.matrix(y[idx]))        
        fit = lm(y ~ x1 + x2, data = testData)
        
        # save all parameters needed for future
        tmp = summary(fit)$coefficients
        predict[guessedMix] = tmp[2, 1]
        predict2[guessedMix] = tmp[3, 1]
        ts[guessedMix]   = tmp[3, 3]
        pval[guessedMix] = tmp[3, 4]

		    # find the candidate contaminator.        
        if (ts[guessedMix] > tsSaved) {
          fitSaved = fit
          tsSaved = ts
          mixSaved = guessedMix
        }
      }
    }
  }
  
  # find the most likely contaminator
  mdx = which.max( ts )

  # print-out
  cat( sprintf( "The possible mixture is %d (%s), with proportion = %.1f%%, t-stat = %.1f, p-value %e\n", 
               mdx, cellNames[mdx], predict2[ mdx ] * 100, ts[mdx], pval[mdx] ) )
               
  return( list( fitSaved, tsSaved, mixSaved )  )
}
