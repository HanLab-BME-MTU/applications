function permXcorrMean = timePermDistXcorrMean(Actmap1, Actmap2, numPerm, lagMax)
% timePermDistXcorrMean Compute time-wise permuation distribution of a xcorr
% curve using parfor by using nanXcorrMaps.m function.
%
% Jungsik Noh, 2016/10/21

% permute time points
tlength = size(Actmap1, 2);
permXcorrMean = nan(numPerm, 2*lagMax+1);
%
parfor i = 1:numPerm;
    permInd = randsample(tlength, tlength);
    
    chan1map = Actmap1(:, permInd);
    chan2map = Actmap2;
        
    xcorrMat_i  = nanXcorrMaps(chan1map, chan2map, 'lagMax', lagMax);
    xcorrMean_i = nanmean(xcorrMat_i, 1);
    
    permXcorrMean(i, :) = xcorrMean_i;
end

