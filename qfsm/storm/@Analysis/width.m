function meanSigmaDist = width(obj,display)

% Compute the width of all the clusters
fun = @(a,b) sqrt(1/numel(a)*sum(sum((b./obj.data.error(a,:)).^2)));
sigmaClusters = cellfun(fun,obj.data.clusters,obj.data.modelRes); % RMS of the distances

sigmaClusters = sigmaClusters/sqrt(2)*9;

if numel(sigmaClusters) > 1
    % Compute the weighted histogram
    % sigmaClusterHistogramWeigths = obj.data.modelLength;
    sigmaClusterHistogramWeigths = obj.data.clusterSize;
    nBins = ceil(sqrt(numel(sigmaClusters)));

    %     binStart = min(sigmaClusters);
    %     binEnd = max(sigmaClusters);
    binStart = 0;
    binEnd = 30;
    
    %     nBins = nBins*4; % 1,2,4,8
    %     binSize = (binEnd-binStart)/nBins
    
    binSize = 0.25
    nBins = floor((binEnd-binStart)/binSize);
    
    histEdges = binStart:binSize:binEnd;
    histValues = zeros(1,nBins);
    for b=1:nBins-1
        binIdx = sigmaClusters>=histEdges(b) & sigmaClusters<histEdges(b+1);
        histValues(1,b) = sum(sigmaClusterHistogramWeigths(binIdx));
    end
    b = nBins;
    binIdx = sigmaClusters>=histEdges(b) & sigmaClusters<=histEdges(b+1);
    histValues(1,b) = sum(sigmaClusterHistogramWeigths(binIdx));
    binCenters = binStart+binSize/2:binSize:binEnd-binSize/2;
    
    % Plot the histogram
    if nargin > 1
        if strcmp(display,'on')
            bar(binCenters,histValues,1);
            xlim([0,3]);
            ylim([0,4000]);
            % xlabel('RMS of the normalized distances');
            xlabel('Standard deviation of the distance components');
            ylabel('Number of points');
        end
    else
        bar(binCenters,histValues,1);
        % xlim([0,3]);
        % ylim([0,4000]);
        % xlabel('RMS of the normalized distances');
        xlabel('Standard deviation of the distance components');
        ylabel('Number of points');
    end
    
else
    fprintf('Analysis: Only one cluster: Normalized RMS distance: %.2f\n',sigmaClusters);
end

meanSigmaDist = sum(sigmaClusters.*obj.data.clusterSize)/sum(obj.data.clusterSize);

end

