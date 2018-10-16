function testStats = CCLE_Identification( ccle_freq, ccle_dp, test_freq, test_dp, targetID, removeCommonFlag )
% Report cell-line identification.
%
%   report = CCLE_Identification( freq, dp, freq(:, 100), dp(:, 100) );
%
%     ccle_freq:    SNP frequency of all cell-lines in database
%     ccle_dp:     	SNP depth of coverage.
%     test_freq:    SNP frequency of test cell-line
%     test_dp:      SNP depth of coverage.
%	  targetID:			optional, only if you have a target cell-line for stats only
%	  removeCommonFlag: if provided, the common SNPs will be removed. See the code. 
%
%   for freq and dp, these are structures with
%       freq.names:    all names of samples. G28824.HEC-265.3, G28573.OV56.1, etc
%       freq.Position: SNP position (1:19505549, 1:19505583) the first number 
%                        is the chr, and second number genomic position).
% .     freq.value:    matrix of all frequency (or dp). 
%

%   Yidong Chen, Jan 2018, UTHSCSA/GCCRI

if( nargin < 6 ), removeCommonFlag = 0; end
if( nargin < 7 ), celllineInfoFileName = 'CCLE_sample_info_file_2012-10-18.txt'; end
    
numSamples = length( ccle_freq.names );

% make sure both ccle_freq and test_freq match.
[c, ia, ib] = intersect( ccle_freq.Position, test_freq.Position );

ccle_freq.Position  = ccle_freq.Position( ia );
ccle_freq.value     = ccle_freq.value( ia, : );
ccle_dp.Position    = ccle_dp.Position( ia );
ccle_dp.value       = ccle_dp.value( ia, : );

test_freq.Position  = test_freq.Position( ib );
test_freq.value     = test_freq.value( ib, : );
test_dp.Position    = test_dp.Position( ib );
test_dp.value       = test_dp.value( ib, : );

if( nargin < 5 )
    fprintf( 1, 'Total of %d mutations selected. \n', length(ccle_freq.Position) );
    testStats = CCLE_CCLE_Identification_allSamples( ccle_freq, ccle_dp, test_freq, test_dp );
elseif( nargin < 6 )
    % calculate corelation coeffcient for targetID in CCLE sample (again the test
    % sample)
    dpThr = 10;
    i = targetID;
    tmp1 = ccle_freq.value(:,i) > 0;
    tmp2 = ccle_dp.value(:,i) >= dpThr;
    snpCounts(i) = sum( tmp1 .* tmp2 );
    idx = find( ( (test_freq.value > 0) | tmp1 ) & ( test_dp.value >= dpThr ) & tmp2 );
    jcount = length(idx);
    c = corrcoef( ccle_freq.value(idx,i), test_freq.value(idx) );
    rho = c(1, 2);
    
    p = normcdf( rho, 0.470, 0.047, 'upper' );
    if( p > 0.001 )
        fprintf( 1, 'match cell-line = %s, with rho = %.2f, p = %.3f\n', ccle_freq.names{targetID}, rho, p );
    else
        fprintf( 1, 'match cell-line = %s, with rho = %.2f, p = %.3e\n', ccle_freq.names{targetID}, rho, p );
    end
    testStats.name = ccle_freq.names{targetID};
    testStats.p = p;
    testStats.rho = rho;
    testStats.id = targetID;
    testStats.numMuts = jcount;
else
    % remove common SNP first. We define comoon SNPs are those FREQ
    % is less than 10% (either difference between min/max. 
    maxFreq = max( ccle_freq.value, [], 2 );
    minFreq = min( ccle_freq.value, [], 2 );
    minDP = min( ccle_dp.value, [], 2 );
    idx = find( ( maxFreq - minFreq ) > 10 & minDP > 20 );

    ccle_freq.Position  = ccle_freq.Position( idx );
    ccle_freq.value     = ccle_freq.value( idx, : );
    ccle_dp.Position    = ccle_dp.Position( idx );
    ccle_dp.value       = ccle_dp.value( idx, : );
    
    test_freq.Position  = test_freq.Position( idx );
    test_freq.value     = test_freq.value( idx, : );
    test_dp.Position    = test_dp.Position( idx );
    test_dp.value       = test_dp.value( idx, : );
    fprintf( 1, 'Total of %d mutations selected. \n', length(idx) );

    testStats = CCLE_CCLE_Identification_allSamples( ccle_freq, ccle_dp, test_freq, test_dp );
end

return

function testStats = CCLE_CCLE_Identification_allSamples( ccle_freq, ccle_dp, test_freq, test_dp )

numSamples = length( ccle_freq.names );

% calculate corelation coeffcient for all ccle samples (again the test
% sample)
rho = ones( numSamples, 1 );
jcount = zeros( numSamples, 1 );
dpThr = 10;

tmp1 = ccle_freq.value > 0;
tmp2 = ccle_dp.value >= dpThr;
t=cputime;
for i = 1:numSamples
    snpCounts(i) = sum( tmp1(:,i) .* tmp2(:,i) );
    idx = find( ( (test_freq.value > 0) | tmp1(:,i) ) & ( test_dp.value >= dpThr ) & tmp2(:,i) );
    jcount(i) = length(idx);
    if( jcount(i) <= 1 )
        rho(i) = 0;
    else
        c = corrcoef( ccle_freq.value(idx,i), test_freq.value(idx) );
        rho( i ) = c(1, 2);
    end
end

% output each cell-line and it best and second to best cell-line name and
% correlation coefficient/number of matching mutations.
[srho, sdx] = sort( rho, 'descend' );

% Within C group (mu, sigma, p0.001) = 0.450, 0.029, 0.539
% Within G group (mu, sigma, p0.001) = 0.470, 0.047, 0.616
% Between C&G groups (mu, sigma, p0.001) = 0.275, 0.042, 0.407
%
% this may have to change.
p = 0;
k = 0;
while( (p < 0.001 && k < numSamples) || k < 2 )
    k = k+1;
    p = normcdf( srho(k), 0.470, 0.047, 'upper' );
    if( k == 1 ) txtStr = '1st'; end
    if( k == 2 ) txtStr = '2nd'; end
    if( k == 3 ) txtStr = '3rd'; end
    if( k > 3 ), txtStr = [num2str(i) 'th']; end
    if( p > 0.001 )
        fprintf( 1, '%s match cell-line = %s, with rho = %.2f, p = %.3f, (%d)\n', txtStr, ccle_freq.names{sdx(k)}, srho(k), p, sdx(k) );
    else
        fprintf( 1, '%s match cell-line = %s, with rho = %.2f, p = %.3e, (%d)\n', txtStr, ccle_freq.names{sdx(k)}, srho(k), p, sdx(k) );
    end
    testStats(k).name = ccle_freq.names{sdx(k)};
    testStats(k).p = p;
    testStats(k).rho = srho(k);
    testStats(k).id = sdx(k);
    testStats(k).numMuts = jcount(sdx(k));
end
if( length( testStats ) == 1 )
    testStats(2) = testStats(1);
end
fprintf( 1, 'Test sample = %s, done in %.2f mins. \n', test_freq.names{1}, (cputime - t )/60 );

return
