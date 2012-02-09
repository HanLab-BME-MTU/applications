function meanSigmaDist = width(obj,display)

% Compute the width of all the clusters
fun = @(a,b) sqrt(1/numel(a)*sum(sum((b./obj.data.error(a,:)).^2)));
sigmaClusters = cellfun(fun,obj.data.clusters,obj.data.modelRes); % RMS of the distances

if numel(sigmaClusters) > 1
    % Compute the weighted histogram
    % sigmaClusterHistogramWeigths = obj.data.modelLength;
    sigmaClusterHistogramWeigths = cellfun(@numel,obj.data.clusters);
    nBins = ceil(sqrt(numel(sigmaClusters)))*4;
    binSize = (max(sigmaClusters)-min(sigmaClusters))/nBins;
    histEdges = min(sigmaClusters):binSize:max(sigmaClusters);
    histValues = zeros(1,nBins);
    for b=1:nBins-1
        binIdx = sigmaClusters>=histEdges(b) & sigmaClusters<histEdges(b+1);
        histValues(1,b) = sum(sigmaClusterHistogramWeigths(binIdx));
    end
    b = nBins;
    binIdx = sigmaClusters>=histEdges(b) & sigmaClusters<=histEdges(b+1);
    histValues(1,b) = sum(sigmaClusterHistogramWeigths(binIdx));
    binCenters = min(sigmaClusters)+binSize/2:binSize:max(sigmaClusters)-binSize/2;
    
    % Plot the histogram
    if nargin > 1
        if strcmp(display,'on')
            bar(binCenters,histValues,1);
            xlim([0,3]);
            ylim([0,4000]);
            xlabel('RMS of the normalized distances');
            ylabel('Number of points');
        end
    else
        bar(binCenters,histValues,1);
        xlim([0,3]);
        ylim([0,4000]);
        xlabel('RMS of the normalized distances');
        ylabel('Number of points');
    end
    
else
    fprintf('Analysis: Only one cluster: Normalized RMS distance: %.2f\n',sigmaClusters);
end

clusterSize = cellfun(@numel,obj.data.clusters);
meanSigmaDist = sum(sigmaClusters.*clusterSize)/sum(clusterSize);

end

