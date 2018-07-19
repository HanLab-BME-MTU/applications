clc
clear
close all

% generate three synthetic gaussian clusters
numPointsPerCluster = 2^10;
ptTrueClusterCenters = [ 1 1 ; -1, -1; 1, -1 ];
clusterStdDev = 0.6;

ptRawData = [];
for i = 1:size(ptTrueClusterCenters,1)
    
    ptCurCluster = repmat(ptTrueClusterCenters(i,:), numPointsPerCluster, 1) + ...
                   clusterStdDev * randn(numPointsPerCluster,2);               
    ptRawData = [ ptRawData; ptCurCluster ];
    
end

% Run mean-shift
kBandwidth = round( 0.3 * size(ptRawData,1) );
tic
profile on;
[clusterInfo,pointToClusterMap] = AdaptiveMeanShiftClustering( ptRawData, kBandwidth, ... 
                                                               'minClusterDistance', 1.8 * clusterStdDev, ... 
                                                               'maxBandwidth', 1.5 * clusterStdDev, ...
                                                               'flagDebug', false, ...
                                                               'kernel', 'gaussian', ...
                                                               'flagUseKDTree', true);
profile off;
profile viewer;
toc

% plot result
figure;
hold on;

    clusterColors = mat2gray( squeeze( label2rgb( 1:numel(clusterInfo) ) ), [0,255] );
    if numel(clusterInfo) == 1
        clusterColors = clusterColors';
    end
    for k = 1:numel(clusterInfo)
        ptCurClusterCenter = mean( ptRawData(clusterInfo(k).ptIdData, :), 1 );
        plot( ptRawData(pointToClusterMap==k, 1), ...
              ptRawData(pointToClusterMap==k, 2), ...
              '.', 'Color', clusterColors(k,:))
        plot(ptCurClusterCenter(1),ptCurClusterCenter(2),'o','MarkerEdgeColor','k','MarkerFaceColor',clusterColors(k,:), 'MarkerSize',10)
    end
    title( sprintf('Mean-shift clustering result - %d clusters were found', numel(clusterInfo)));

hold off;
