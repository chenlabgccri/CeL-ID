function [ccle_id, ccle_name, prcTest, prcMixture, mdlSaved] = CCLE_MixtureEstimate( ccle_freq, ccle_dp, test_data, test_dp, matched_ccleID )


dpThr = 10;
numSamples = length( ccle_freq.names );

predict = zeros( numSamples, 1 );
predict2 = zeros( numSamples, 1 );
ts = zeros( numSamples, 1 );
tsSaved = 0;
for guessedMix = 1:numSamples
    x1 = ccle_freq.value(:,matched_ccleID);
    x2 = ccle_freq.value(:,guessedMix);
    y = test_data.value;
    
    dp1 = ccle_dp.value(:,matched_ccleID);
    dp2 = ccle_dp.value(:,guessedMix);
    idx = find( dp1 >= dpThr & dp2 >= dpThr );
    
    x1 = x1(idx);
    x2 = x2(idx);
    y = y(idx);

    aa = x1 - x2;
    [saa, adx] = sort( aa );
    
    numMuts = length( aa );
    midPoint = fix( numMuts/2 );
    idx = [adx(1:200); adx(midPoint-100:midPoint+100); adx(end-200:end)];
    mdl = fitlm( [x1(idx) x2(idx)], y(idx), 'linear', 'RobustOpts', 'on' );
    
    params = mdl.Coefficients.Estimate;
    tstats = mdl.Coefficients.tStat;
    tpval = mdl.Coefficients.pValue;
    predict( guessedMix ) = params(2);
    predict2( guessedMix ) = params(3);
    ts( guessedMix ) = tstats(2);
    pval( guessedMix ) = tpval(2);
    
    if( tstats(2) > tsSaved )
        mdlSaved = mdl;
        tsSaved = tstats(2);
    end
end

[mts, mdx] = max( ts );
ccle_id = mdx;
ccle_name = ccle_freq.names{ccle_id};
prcTest = predict2( mdx );
prcMixture = predict( mdx );

fprintf( 1, 'The possible mixture is %d (%s), with proportion = %.1f%%, with t-stat = %.1f, p-value %e\n', ...
            ccle_id, ccle_name, prcMixture * 100, ts(mdx), pval(mdx) );
return