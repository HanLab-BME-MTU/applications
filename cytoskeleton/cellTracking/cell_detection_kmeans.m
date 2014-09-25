Mask = MD.processes_{1}.loadChannelOutput(1,1);

[inda,indb] = find(Mask>0);
% nClass = 6;
% Kmeans_label = kmeans([inda indb],nClass,'distance','hamming','start','cluster');
% 
% figure(1); hold off;
% imagescc(Mask); hold on;
% 
% for iClass = 1 : nClass
%     plot(indb(find(Kmeans_label==iClass)),inda(find(Kmeans_label==iClass)),'.','color',color_array(iClass,:));
% end

ptRawData = [inda indb];
% Run mean-shift
bandwidth = 100;
tic
profile on;
[clusterInfo,pointToClusterMap] = MeanShiftClustering(ptRawData, bandwidth, ... 
                                                      'flagDebug', false, ...
                                                      'kernel', 'gaussian', ...
                                                      'flagUseKDTree', true);
profile off;
profile viewer;
toc

% plot result
figure(2);
hold on;

    clusterColors = mat2gray( squeeze( label2rgb( 1:numel(clusterInfo) ) ), [0,255] );
    for k = 1:numel(clusterInfo)
        ptCurClusterCenter = clusterInfo(k).ptClusterCenter;
        plot( ptRawData(pointToClusterMap==k, 1), ...
              ptRawData(pointToClusterMap==k, 2), ...
              '.', 'Color', clusterColors(k,:));
        sigma_x(k) = std(ptRawData(pointToClusterMap==k, 1));
        sigma_y(k) = std(ptRawData(pointToClusterMap==k, 2));
        
        plot(ptCurClusterCenter(1),ptCurClusterCenter(2),'o','MarkerEdgeColor','k','MarkerFaceColor',clusterColors(k,:), 'MarkerSize',10)
        
    end
    title( sprintf('Mean-shift clustering result - %d clusters were found', numel(clusterInfo)));

hold off;



