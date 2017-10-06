clc
clear
close all

% generate three synthetic gaussian clusters
numPointsPerCluster = 2^12;
ptTrueClusterCenters = [ 1 1 ; -1, -1; 1, -1 ];
clusterStdDev = 0.6;

ptRawData = [];
for i = 1:size(ptTrueClusterCenters,1)
    
    ptCurCluster = repmat(ptTrueClusterCenters(i,:), numPointsPerCluster, 1) + ...
                   clusterStdDev * randn(numPointsPerCluster,2);               
    ptRawData = [ ptRawData; ptCurCluster ];
    
end

% Run mean-shift
bandwidth = 0.5;
tic
numRuns = 100;

h = waitbar(0, 'Running multiple runs of mean-shift ...');

for i = 1:numRuns    
    
    [clusterInfo,pointToClusterMap] = MeanShiftClustering(ptRawData, bandwidth, ... 
                                                          'flagDebug', false, ...
                                                          'kernel', 'gaussian', ...
                                                          'flagUseKDTree', true);
                                                      
    waitbar(i/numRuns, h, sprintf( 'Running multiple runs of mean-shift ... %.2f%% done', 100*i/numRuns) );
    
end

fprintf('\nIt took an average of %.2f seconds to cluster %d points\n', toc / numRuns, size(ptRawData,1));

close(h);

% plot result
figure;
hold on;

    clusterColors = mat2gray( squeeze( label2rgb( 1:numel(clusterInfo) ) ), [0,255] );
    for k = 1:numel(clusterInfo)
        ptCurClusterCenter = clusterInfo(k).ptClusterCenter;
        plot( ptRawData(pointToClusterMap==k, 1), ...
              ptRawData(pointToClusterMap==k, 2), ...
              '.', 'Color', clusterColors(k,:))
        plot(ptCurClusterCenter(1),ptCurClusterCenter(2),'o','MarkerEdgeColor','k','MarkerFaceColor',clusterColors(k,:), 'MarkerSize',10)
    end
    title( sprintf('Mean-shift clustering result - %d clusters were found', numel(clusterInfo)));

hold off;
