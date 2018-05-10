function CCLE_Stats( FREQfileName, DPfileName, DNAFREQfileName, DNADPfileName, testID,  celllineInfoFileName )
% Get statistics for cell line through correlation coefficient calculation
%
%   stats = CCLE_Stats( 'CCLE_v2_filtered_FREQ.txt', 'CCLE_v2_filtered_DP.txt' );
%   stats = CCLE_Stats( 'CCLE_v2_filtered_FREQ.txt', 'CCLE_v2_filtered_DP.txt', 4,  'CCLE_sample_info_file_2012-10-18.txt' );
%
%     FREQfileName: file name to SNP frequency
%     DPfileName:   file name to SNP depth of coverage.
%     celllineInfoFileName: (optional) it is a cell-line information, include cell-line names, etc.
%
%   both file are in matrix format: first row is CCLE cell-line name (G28824.HEC-265.3, G28573.OV56.1)
%   first row is SNP position (1:19505549, 1:19505583, the first number is the chr, and second number
%   genomic position). 
%

%   Yidong Chen, Jan 2018, UTHSCSA/GCCRI

if( nargin < 5 ), testID = 33; end
if( nargin < 6 ), celllineInfoFileName = 'CCLE_sample_info_file_2012-10-18.txt'; end
    
data = ReadTextWithHeader( FREQfileName, 3 );
dp = ReadTextWithHeader( DPfileName, 3 );

numSamples = length( data.names );
rho = ones( numSamples, numSamples );
jcount = zeros( numSamples, numSamples );
dpThr = 10;

tmp1 = data.value > 0;
tmp2 = dp.value >= dpThr;
t=cputime;
for i = 1:numSamples
    snpCounts(i) = sum( tmp1(:,i) .* tmp2(:,i) );
    jcount(i,i) = snpCounts(i);
    if( mod( i, 100 ) == 0 ), fprintf( 1, '%d (%.1f).', i, (cputime-t)/60 ); end;
    a1 = tmp1( :, i );
    a2 = tmp2( :, i);
    for j = 1:numSamples
        if( j > i )
            idx = find( ( a1 | tmp1(:,j) ) & a2 & tmp2(:,j) );
            jcount(i, j ) = length(idx);
            jcount(j, i ) = jcount(i, j );
            c = corrcoef( data.value(idx,i), data.value(idx, j) );
            rho( i, j ) = c(1, 2);
            rho( j, i ) = c(1, 2);
        end
    end
end
fprintf( 1, 'done. \n' );

% generate dendrogram
Y = squareform( 1-rho );
Z = linkage( Y, 'average' );
[H, T, perm] = dendrogram( Z, 0, 'label', data.names );
ylabel( '1 - correlation' );
rotateticklabel( gca, 90 );
orient landscape
print( '-dpdf', 'CCLE-dendrogram-correlation.pdf' );

fid = fopen( 'CCLE-dendrogram-correlation-sampleNames.txt', 'w' );
for i = 1:length( perm )
    fprintf( fid, '%s\n', data.names{perm(i)} );
end
fclose( fid );


% output each cell-line and it best and second to best cell-line name and
% correlation coefficient/number of matching mutations. 
fid = fopen( 'CCLE_stats_perCellLines.txt', 'w' );
fprintf( fid, 'order\tcell-line name\tNumber of mutations\t1st match CL\t1st-corr\t1st-snpCount\t2nd match CL\t2nd-corr\t2nd-snpCount\n' );
for i = 1:numSamples
    fprintf( fid, '%d\t%s\t%d\t', i, data.names{i}, jcount(i,i)  );
    
    [srho, sdx] = sort( rho(i,:), 'descend' );
    fprintf( fid, '%s\t%.3f\t%d\t', data.names{sdx(2)}, rho(i, sdx(2)), jcount(i, sdx(2)) );
    fprintf( fid, '%s\t%.3f\t%d\n', data.names{sdx(3)}, rho(i, sdx(3)), jcount(i, sdx(3)) );
    
    if( rho(i, sdx(2)) > 0.8 ), fprintf( 1, '%d, %s\n', i, data.names{i} ); end
end
fclose( fid );


% generate a test figure;
figure;
testID = 33;
subplot( 3, 1, 1 );
[srho, sdx] = sort( rho(testID,:), 'descend' );
rr = rho( :, testID );
rr( testID ) = [];
hist( rr, 50 );
xlabel( ['correlation with ' data.names{testID}], 'interpret', 'none' );
text( rho(testID, sdx(2)), 10, [data.names{sdx(2)} ' (' num2str( rho(testID, sdx(2)), '%.2f') ')'], ...
      'horizontalalignment', 'left', 'verticalalignment', 'middle', 'rotation', 90, 'interpret', 'none' );
ylabel( 'Number of samples' );

subplot( 3, 1, 2 );
jc = jcount( :, testID );
hist( jc/1000, 50 );
axis( [0 max(jc/1000) 0 Inf] );
xlabel( ['# of common mutations with ' data.names{testID} ' (1,000)'], 'interpret', 'none' );
text( jcount(testID, sdx(2))/1000, 10, [data.names{sdx(2)} ' (' num2str( jcount(testID, sdx(2))) ')'], ...
      'horizontalalignment', 'left', 'verticalalignment', 'middle', 'rotation', 90, 'interpret', 'none' );
ylabel( 'Number of samples' );

subplot( 3, 2, 5 );
idx = find( ( tmp1(:,testID) | tmp1(:,sdx(2)) ) & tmp2(:,testID) & tmp2(:,sdx(2)) );
plot( data.value( idx, testID ), data.value( idx, sdx(2) ), '.' ); 
xlabel( [data.names{testID} ', FREQ (%)'], 'interpret', 'none' );
ylabel( [data.names{sdx(2)} ', FREQ (%)'], 'interpret', 'none' );
title( ['rho = ' num2str( rho(testID, sdx(2)) )] );
refline( 1, 0 );

subplot( 3, 2, 6 );
idx = find( ( tmp1(:,testID) | tmp1(:,sdx(3)) ) & tmp2(:,testID) & tmp2(:,sdx(3)) );
plot( data.value( idx, testID ), data.value( idx, sdx(3) ), '.' ); 
xlabel( [data.names{testID} ', FREQ (%)'], 'interpret', 'none' );
ylabel( [data.names{sdx(3)} ', FREQ (%)'], 'interpret', 'none' );
title( ['rho = ' num2str( rho(testID, sdx(3)) )] );
refline( 1, 0 );

orient portrait
orient tall
print( '-dpdf', ['CCLE-stats-' data.names{testID} '.pdf'] );


% we model correlation coefficient as normal distribution. Howevever, in
% order to avoid outliers, we use the truncated normal with r < 0.8.
% Therefore the truncated normal distribution is as follows,
norm_trunc =@(x, mu, sigma) (normpdf(x , mu, sigma)./normcdf(0.8, mu, sigma) .* heaviside(0.8-x));
% The truncated pdf should be normilized. since we truncate right side, then 
% divide by normcdf( 0.8, mu, sigma) 

% find distribution of correlation coefficient.
siteFlag = zeros( length( data.names ), 1 );
for i = 1:length( siteFlag )
    if( data.names{i}(1) == 'C' ), siteFlag(i) = 1; end
end

% keep all correlation within site. 
sf1 = repmat( siteFlag, 1, length( data.names ) );
sf2 = repmat( siteFlag', length( data.names ), 1 );
rr = rho .* sf1 .* sf2;
rr = (-2)*tril(ones(size(rr)), -1) + triu( rr );
rr = rr(:);
rr = rr( rr < 1 & rr > 0 );
[n1, x1] = hist( rr, 0:0.02:1 );

%find the maximum likelihood estimates using mean and std as an initial guess 
rr = rr( rr < 0.8 );
[phat1, phat1_ci] = mle(rr, 'pdf', norm_trunc, 'start', [mean(rr), std(rr)] );
fprintf( 1, 'Within WES samples (mu, sigma, p0.001) = %.3f, %.3f, %.3f\n', phat1(1), phat1(2), norminv( 0.999, phat1(1), phat1(2) ) );


% now we know sample name start with C are DNA samples, let's load it.
%dna_data = ReadTextWithHeader( DNAFREQfileName, 3 );
%dna_dp = ReadTextWithHeader( DNADPfileName, 3 );

sf1 = repmat( ~siteFlag, 1, 1260 );
sf2 = repmat( ~siteFlag', 1260, 1 );
rr = rho .* sf1 .* sf2;
rr = (-2)*tril(ones(size(rr)), -1) + triu( rr );
rr = rr(:);
rr = rr( rr < 1 & rr > 0 );
[n2, x2] = hist( rr, 0:0.02:1 );

%find the maximum likelihood estimates using mean and std as an initial guess 
rr = rr( rr < 0.8 );
[phat2, phat2_ci] = mle(rr, 'pdf', norm_trunc, 'start', [mean(rr), std(rr)] );
fprintf( 1, 'Within RNA-seq samples (mu, sigma, p0.001, p0.0001) = %.3f, %.3f, %.3f, %.3f\n', ...
        phat2(1), phat2(2), norminv( 0.999, phat2(1), phat2(2) ), norminv( 0.999999, phat2(1), phat2(2) ) );

sf1 = repmat( siteFlag, 1, length( data.names ) );
sf2 = repmat( ~siteFlag', length( data.names ), 1 );
rr = rho .* sf1 .* sf2;
rr = (-2)*tril(ones(size(rr)), -1) + triu( rr );
rr = rr(:);
rr = rr( rr < 1 & rr > 0 );
[n3, x3] = hist( rr, 0:0.02:1 );

%find the maximum likelihood estimates using mean and std as an initial guess 
rr = rr( rr < 0.8 );
[phat1, phat1_ci] = mle(rr, 'pdf', norm_trunc, 'start', [mean(rr), std(rr)] );
fprintf( 1, 'Between WES & RNAseq (mu, sigma, p0.001) = %.3f, %.3f, %.3f\n', phat1(1), phat1(2), norminv( 0.999, phat1(1), phat1(2) ) );

figure;
subplot( 2, 1, 2 );
hdl1 = bar( x1, n1/1000 );
hold on;
hdl2 = bar( x2, n2/1000 );
hdl3 = bar( x3, n3/1000 );
hdl4 = plot( x2, normpdf( x2, phat2(1), phat2(2) )*(sum(n2)*(x2(2)-x2(1)))/1000, 'b' );
set( hdl4, 'linewidth', 2 )

xlabel( 'correlation' );
ylabel( '# of samples (x1000)' );
set( hdl2, 'FaceAlpha', 0.5 );
set( hdl3, 'FaceAlpha', 0.5, 'FaceColor', 'r' );
legend( [hdl1 hdl2 hdl3 hdl4], [{'within WES'} {'within RNA'} {'RNAseq vs WES'} {'MLE of corr of RNA'}] );
axis( [0 1 0 ceil(max( [n1 n2 n3]/1000 )/10)*10] );
orient tall
print( '-dpdf', 'CCLE-correlation-distributionExplained.pdf' );
close;

fid = fopen( 'CCLE_clusters_0.8.txt', 'w' );
T = cluster( Z, 'cutoff', 0.8 );
k = 0;
for i = 1:max( T )
    if( length( find( T == i ) ) > 1 )
        dx = find( T == i );
        k = k+1;
        fprintf( fid, 'k = %d, Cluster = %d, samples = ', k, i );
        for j = 1:length(dx)
            fprintf( fid, '{%d, %s},', dx(j), data.names{dx(j)} );
        end
        fprintf( fid, '%.3f, %d', rho(dx(1), dx(2)), jcount(dx(1), dx(2)) );
        fprintf( fid, '\n' );
    end
end
fclose( fid );

% generate correlation and mutation count figure.
figure;
r = rho(:);
jcc = zeros( numSamples, numSamples );
for i = 1:numSamples
    for j = 1:numSamples
        jcc(i,j) = jcount(i,j) * 100 / min( jcount(i,i), jcount(j, j) );
        if( jcc(i, j) > 80 && jcc(i, j) < 100 && j > i )
            fprintf( 1, '%s\t%s\t%.3f\t%d\t%d\n', data.names{i}, data.names{j}, rho(i, j), i, j );
        end
    end
end
jcc = jcc(:);

plot( r, jcc, '.' );
xlabel( 'Correlation' );
ylabel( 'log2( overlap SNP counts )' );
set( gcf, 'papersize', [10 10] );
orient tall
print( '-dpdf',  'CCLE_correlation_mutCount_Figure.pdf' );
close

return