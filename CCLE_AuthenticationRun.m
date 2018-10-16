FREQfileName = 'CCLE_v2_filtered_FREQ.txt';
DPfileName = 'CCLE_v2_filtered_DP.txt';

% this program may take 1 hour to run. It generates all stats needed to
% characterize CCLE mutation info.
CCLE_Stats( FREQfileName, DPfileName, 4 );

% Cell line identification via correlation coefficient.
FREQfileName = 'RNASeq/CCLE_RNA_recode_FREQ.txt';
DPfileName = 'RNASeq/CCLE_RNA_recode_DP.txt';
data = ReadTextWithHeader( FREQfileName, 3 );
dp = ReadTextWithHeader( DPfileName, 3 );

clear accMu accSigma accMu2 accSigma2
% generate test data from CCLE set.
testSampleIDs = [33 278 691 636 109 140];
clear acc acc2;
for i = 1:length( testSampleIDs )
    test_data.names = data.names( testSampleIDs(i) );
    test_data.Position = data.Position;
    test_data.value = data.value( :, testSampleIDs(i) );
    test_dp = test_data;
    test_dp.value = dp.value(:, testSampleIDs(i));
    stats = CCLE_Identification( data, dp, test_data, test_dp );  % takes about 1 min.
    acc(i) = stats(2).rho;
    acc2(i) = stats(3).rho;
end
accMu(1) = mean( acc );
accSigma(1) = std( acc );
accMu2(1) = mean( acc2 );
accSigma2(1) = std( acc2 );

% similar, but using costmic data.
FREQfileName = '/Users/cheny8/Dropbox/StudentFolderShared/CCLE_notshared/RNASeq/Cosmic70_RNA_matched_FREQ.txt';
DPfileName = '/Users/cheny8/Dropbox/StudentFolderShared/CCLE_notshared/RNASeq/Cosmic70_RNA_matched_DP.txt';
data = ReadTextWithHeader( FREQfileName, 3 );
dp = ReadTextWithHeader( DPfileName, 3 );

% generate test data from CCLE set.
testSampleIDs = [33 278 691 636 109 140];
clear acc acc2;
for i = 1:length( testSampleIDs )
    test_data.names = data.names( testSampleIDs(i) );
    test_data.Position = data.Position;
    test_data.value = data.value( :, testSampleIDs(i) );
    test_dp = test_data;
    test_dp.value = dp.value(:, testSampleIDs(i));
    stats = CCLE_Identification( data, dp, test_data, test_dp );  % takes about 1 min.
    acc(i) = stats(2).rho;
    acc2(i) = stats(3).rho;
end
accMu(2) = mean( acc );
accSigma(2) = std( acc );
accMu2(2) = mean( acc2 );
accSigma2(2) = std( acc2 );

% similar, but using costmic 83.
FREQfileName = '/Users/cheny8/Dropbox/StudentFolderShared/CCLE_notshared/RNASeq/Cosmic83_RNA_matched_FREQ.txt';
DPfileName = '/Users/cheny8/Dropbox/StudentFolderShared/CCLE_notshared/RNASeq/Cosmic83_RNA_matched_DP.txt';
data = ReadTextWithHeader( FREQfileName, 3 );
dp = ReadTextWithHeader( DPfileName, 3 );

% generate test data from CCLE set.
testSampleIDs = [33 278 691 636 109 140];
clear acc acc2;
for i = 1:length( testSampleIDs )
    test_data.names = data.names( testSampleIDs(i) );
    test_data.Position = data.Position;
    test_data.value = data.value( :, testSampleIDs(i) );
    test_dp = test_data;
    test_dp.value = dp.value(:, testSampleIDs(i));
    stats = CCLE_Identification( data, dp, test_data, test_dp );  % takes about 1 min.
    acc(i) = stats(2).rho;
    acc2(i) = stats(3).rho;    
end
accMu(3) = mean( acc );
accSigma(3) = std( acc );
accMu2(3) = mean( acc2 );
accSigma2(3) = std( acc2 );


% Use WES data.
FREQfileName = '/Users/cheny8/Dropbox/StudentFolderShared/CCLE_notshared/RNASeq/CCLE_RNA_recode_FREQ.txt';
DPfileName = '/Users/cheny8/Dropbox/StudentFolderShared/CCLE_notshared/RNASeq/CCLE_RNA_recode_DP.txt';
data = ReadTextWithHeader( FREQfileName, 3 );
dp = ReadTextWithHeader( DPfileName, 3 );

% generate test data from CCLE set.
FREQfileName = '/Users/cheny8/Dropbox/StudentFolderShared/CCLE_notshared/WESeq/CCLE_WES_recode_FREQ.txt';
DPfileName = '/Users/cheny8/Dropbox/StudentFolderShared/CCLE_notshared/WESeq/CCLE_WES_recode_DP.txt';
tdata = ReadTextWithHeader( FREQfileName, 3 );
tdp = ReadTextWithHeader( DPfileName, 3 );

clear acc acc2;
for i = 1:length( tdata.names )
    test_data.names = tdata.names(i);
    test_data.Position = tdata.Position;
    test_data.value = tdata.value( :, i );
    test_dp = test_data;
    test_dp.value = tdp.value(:, i);
    stats = CCLE_Identification( data, dp, test_data, test_dp );  % takes about 1 min.
    acc(i) = stats(1).rho;
    if( stats(2).rho > .8 )
        acc2(i) = stats(3).rho;
    else
        acc2(i) = stats(2).rho;
    end
end
accMu(4) = mean( acc );
accSigma(4) = std( acc );
accMu2(4) = mean( acc2 );
accSigma2(4) = std( acc2 );

figure;
bar(1:4,accMu2)
hold on
errorbar(1:4,accMu2,accSigma2,'.')

% Scramble...
FREQfileName = '/Users/cheny8/Dropbox/StudentFolderShared/CCLE_notshared/RNASeq/CCLE_RNA_recode_FREQ.txt';
DPfileName = '/Users/cheny8/Dropbox/StudentFolderShared/CCLE_notshared/RNASeq/CCLE_RNA_recode_DP.txt';
data = ReadTextWithHeader( FREQfileName, 3 );
dp = ReadTextWithHeader( DPfileName, 3 );

% generate test data from CCLE set.
testSampleIDs = [33 278 691 636 109 140];
clear acc;
for i = 1:length( testSampleIDs )
    test_data.names = data.names( testSampleIDs(i) );
    test_data.Position = data.Position( randperm(length(data.Position)) );
    test_data.value = data.value( :, testSampleIDs(i) );
    test_dp = test_data;
    test_dp.value = dp.value(:, testSampleIDs(i));
    stats = CCLE_Identification( data, dp, test_data, test_dp );  % takes about 1 min.
    acc(i) = stats(1).rho;
end
accMu(5) = mean( acc );
accSigma(5) = std( acc );

figure;
bar(1:5,accMu);
hold on
errorbar(1:5,accMu,accSigma,'.');
set(gca, 'xticklabel', [{'RNAseq'} {'RNAseq-Cosmic70'} {'RNAseq-Cosmic70'} {'RNAseq-WES'} {'randomPerm'}] );


% robustness test. 
% Cell line identification via correlation coefficient.
FREQfileName = '/Users/cheny8/Dropbox/StudentFolderShared/CCLE_notshared/RNASeq/CCLE_RNA_recode_FREQ.txt';
DPfileName = '/Users/cheny8/Dropbox/StudentFolderShared/CCLE_notshared/RNASeq/CCLE_RNA_recode_DP.txt';
data = ReadTextWithHeader( FREQfileName, 3 );
dp = ReadTextWithHeader( DPfileName, 3 );

% generate test data from CCLE set.
testSampleIDs =   [ 33 278 691 636 109 140];
targetSampleIDs = [686 933 932 869 129 717];
col = 'rgbcmk';
for i = 1:length( testSampleIDs )
    clear acc;
    nLevels = [1 5 10 20 30 50];
    test_data.names     = data.names( testSampleIDs(i) );
    test_data.Position  = data.Position;
    test_data.value     = data.value(:, testSampleIDs(i));
    test_dp.names       = dp.names( testSampleIDs(i) );
    test_dp.Position    = dp.Position;
    test_dp.value       = dp.value(:, testSampleIDs(i));

    target_data = data;
    target_data.value = target_data.value(:, targetSampleIDs(i));
    target_data.names = target_data.names(targetSampleIDs(i));
    target_dp = dp;
    target_dp.names = target_dp.names(targetSampleIDs(i));
    target_dp.value = target_dp.value(:, targetSampleIDs(i));
    for j = 1:length( nLevels )
        clear acc;
        for jj = 1:10
            test_data.value = normrnd( data.value( :, testSampleIDs(i) ), nLevels(j) );
            test_data.value( test_data.value < 0 ) = 0;
            test_data.value( test_data.value > 100 ) = 100;
            stats = CCLE_Identification( target_data, target_dp, test_data, test_dp, 1 );  % takes about 1 min.
            acc(jj) = stats.rho;
        end
        accMu(j) = mean( acc );
        accSigma(j) = std( acc );
    end
    plot( nLevels, accMu, col(i) )
    if( i == 1 ), hold on; end
    errorbar( nLevels, accMu, accSigma, '.' );
end
ylabel( 'Correlation' );
xlabel( 'noise level' );
print -dpdf Robustness_6pairedSamples.pdf

% Test MCF7 samples against CCLE RNAseq samples
FREQfileName = '/Users/cheny8/Dropbox/StudentFolderShared/CCLE_notshared/RNASeq/CCLE_RNA_recode_FREQ.txt';
DPfileName = '/Users/cheny8/Dropbox/StudentFolderShared/CCLE_notshared/RNASeq/CCLE_RNA_recode_DP.txt';
data = ReadTextWithHeader( FREQfileName, 3 );
dp = ReadTextWithHeader( DPfileName, 3 );

% Get samples from GSE86316
FREQfileName = '/Users/cheny8/Dropbox/StudentFolderShared/CCLE_notshared/testData/GSE86316_FREQ_nochr.txt';
DPfileName = '/Users/cheny8/Dropbox/StudentFolderShared/CCLE_notshared/testData/GSE86316_DP_nochr.txt';
tdata = ReadTextWithHeader( FREQfileName, 3 );
tdp = ReadTextWithHeader( DPfileName, 3 );

clear acc;
clear acc2;
clear acc3;
clear ct;
for i = 1:length( tdata.names )
    test_data.names = tdata.names(i);
    test_data.Position = tdata.Position;
    test_data.value = tdata.value( :, i );
    test_dp = test_data;
    test_dp.value = tdp.value(:, i);
    stats = CCLE_Identification( data, dp, test_data, test_dp );  % takes about 1 min.
    acc(i) = stats(1).rho;
    acc2(i) = stats(2).rho;
    acc3(i) = stats(3).rho;
    ct(i) = stats(1).numMuts;
    
    candidateIDs = [904 486];
    ccle_data.names = data.names(candidateIDs);
    ccle_data.Position = data.Position;
    ccle_data.value = data.value( :, candidateIDs );
    ccle_dp = ccle_data;
    ccle_dp.value = dp.value(:, candidateIDs);
    stats = CCLE_Identification( ccle_data, ccle_dp, test_data, test_dp, [], 1 );  % takes about 1 min.    
end
figure;
bar( [acc; acc2; acc3]' );
print -dpdf MCF7_testResults.pdf


% Test Downsampling against CCLE RNAseq samples
FREQfileName = '/Users/cheny8/Dropbox/StudentFolderShared/CCLE_notshared/RNASeq/CCLE_RNA_recode_FREQ.txt';
DPfileName = '/Users/cheny8/Dropbox/StudentFolderShared/CCLE_notshared/RNASeq/CCLE_RNA_recode_DP.txt';
data = ReadTextWithHeader( FREQfileName, 3 );
dp = ReadTextWithHeader( DPfileName, 3 );

% Get samples from 11 downsampling case.
FREQfileName = '/Users/cheny8/.DellEMS/Vols/disk3s2/Dropbox/StudentPostDocFolderShared/CCLE_notshared/testData/11sample_subsampling_FREQ.txt';
DPfileName   = '/Users/cheny8/.DellEMS/Vols/disk3s2/Dropbox/StudentPostDocFolderShared/CCLE_notshared/testData/11sample_subsampling_DP.txt';
tdata = ReadTextWithHeader( FREQfileName, 3 );
tdp = ReadTextWithHeader( DPfileName, 3 );


FREQfileName = '/Users/cheny8/Dropbox/StudentPostDocFolderShared/CCLE_notshared/testData/GSE101966_FREQ.txt';
DPfileName = '/Users/cheny8/Dropbox/StudentPostDocFolderShared/CCLE_notshared/testData/GSE101966_DP.txt';
tdata = ReadTextWithHeader( FREQfileName, 3 );
tdp = ReadTextWithHeader( DPfileName, 3 );

clear acc;
clear acc2;
clear acc3;
clear ct;
for i = 1:length( tdata.names )
    test_data.names = tdata.names(i);
    test_data.Position = tdata.Position;
    test_data.value = tdata.value( :, i );
    test_dp = test_data;
    test_dp.value = tdp.value(:, i);
    stats = CCLE_Identification( data, dp, test_data, test_dp );  % takes about 1 min.
    acc(i) = stats(1).rho;
    acc2(i) = stats(2).rho;
    ct(i) = stats(1).numMuts;
end

accM = zeros( 11, 6 );
accM2 = zeros( 11, 6 );
ctM = zeros( 11, 6 );
for i = 1:11
    accM(i, 1:6)  =  acc(((i-1)*6+1):i*6);
    accM2(i, 1:6) = acc2(((i-1)*6+1):i*6);
    ctM(i, 1:6)   =   ct(((i-1)*6+1):i*6);
end
ct = [30689 26611 31576 27271 26707 31428 29837 7687 9105];
ac = nanmean( accM(2:10, 1:6) );
ac2 = nanmean( accM2(2:10, 1:6) );
cm = nanmean( ctM(2:10, 1:6) ./ repmat( ct', 1, 6 ) );
figure;
yyaxis left
hdl1 = plot( 1:6, ac, 'c' );
hold on;
errorbar(1:6,ac,nanstd(accM(2:10, 1:6)),'.');
hdl2 = plot( 1:6, ac2, 'b' );
errorbar(1:6,ac2,nanstd(accM2(2:10, 1:6)),'.');
ylabel( 'Correlation' );
yyaxis right
hdl3 = plot( 1:6, cm*100, 'r' );
ylabel( 'Percent of mutations utilized' );


% correlation to expression. 
CCLE_expressionName = '/Users/cheny8/Dropbox/StudentPostDocFolderShared/CCLE_notshared/GeneExpression/CCLE_RNAseq_rpkm.txt';
test_expressionName = '/Users/cheny8/Dropbox/StudentPostDocFolderShared/CCLE_notshared/testData/GSE101966_RSEM_genes_FPKM.txt';
ccleExprData = ReadTextWithHeader( CCLE_expressionName, 3 );
tExprData = ReadTextWithHeader( test_expressionName, 3 );

sampleMappingName = '/Users/cheny8/Dropbox/StudentPostDocFolderShared/CCLE_notshared/GeneExpression/SampleNameMapping.txt';
sList = ReadTextWithHeader( sampleMappingName, 1 );
dx = strmatch( 'G20490.HCT_116.2', sList.CCLE_bamName );

[c, ia, ib] = intersect( ccleExprData.Description, tExprData.ProbeIDs );

a = ccleExprData.value(ia, dx);
b = tExprData.value(ib, 1);
idx = find( a > 0.1 & b > 0.1 );
co = corrcoef( log2(a(idx)), log2( b(idx)) );

for j = 1:8
    b = tExprData.value(ib, j);
    for i = 1:length( ccleExprData.names )
        a = ccleExprData.value(ia, i);
        idx = find( a > 1 & b > 1 );
        [co, p] = corr( log2(a(idx)), log2( b(idx)) );
        acc(i) = co;
        acc2(i) = length( idx );
        accp(i) = p;
    end
    [sacc, sdx] = sort( acc, 'descend' );
    fprintf( 1, '%s\t%s\t%.3f\t%d\n', tExprData.names{j}, ccleExprData.names{sdx(1)}, sacc(1), acc2(sdx(1)) );
    fprintf( 1, '%s\t%s\t%.3f\t%d\n', tExprData.names{j}, ccleExprData.names{sdx(2)}, sacc(2), acc2(sdx(2)) );
end


a = ccleExprData.value(ia, dx);
for i = 1:length( tExprData.names )
    b = tExprData.value(ib, i);
    idx = find( a > 0.1 & b > 0.1 );
    co = corrcoef( log2(a(idx)), log2( b(idx)) );
    acc2(i) = co(1,2);
end

% generate mixture data. 
testSampleID = 33;
mixtureID = 10;

p1 = normrnd(0.85, 0.05, length( data.Position ), 1 );
p2 = normrnd(0.15, 0.02, length( data.Position ), 1 );
mixFreq = data.value( :, testSampleID ) .* p1  + data.value(:, mixtureID) .* p2;
test_data.names{1} = 'mix(0.7*33+0.3*10)';
test_data.Position = data.Position;

test_dp.value = fix( dp.value( :, testSampleID ) .* p1  + dp.value(:, mixtureID) .* p2 );
test_dp.value( test_data.value < 0 ) = 0;  % nothing below 0, but could be as big as needed.
test_dp.names{1} = 'mix(0.7*33+0.3*10)';
test_dp.Position = dp.Position;

noiseSigma = [0.01 0.2 0.5 1 2 5 10 15 20];
clear prcM prcT names;
for i = 1:length( noiseSigma )
    test_data.value = normrnd( mixFreq, noiseSigma(i) );
    test_data.value( test_data.value < 15 ) = 0;        % the data shall be between [0 100].
    test_data.value( test_data.value > 100 ) = 100;
    
    % using linear model to estimate the proportion (it takes about 1-2 min).
    [ccle_id, ccle_name, prcTest, prcMixture, mdl] = CCLE_MixtureEstimate( data, dp, test_data, test_dp, testSampleID );
    prcM(i) = prcMixture;
    prcT(i) = mdl.Coefficients.tStat(2);
    names{i} = ccle_name;    
    % it shall generate something like:
    %   The possible mixture is 10 (G27482.SIMA.2), with proportion = 30.0%, with t-stat = 74.5, p-value 0.000000e+00
end


prcM30 = prcM;
prcT30 = prcT;


figure;
yyaxis left
hdl1 = plot( noiseSigma, 70-prcM30*100, '-+' );
hold on;
hdl2 = plot( noiseSigma, 85-prcM*100, '-x' );
ylabel( 'Proportion Estimate Error' );
xlabel( 'Noise level (sigma)' );
yyaxis right
hdl3 = plot( noiseSigma, prcT30, '-+' );
hold on;
hdl4 = plot( noiseSigma, prcT, '-x' );
ylabel( 't-stat of linear model fit' );

