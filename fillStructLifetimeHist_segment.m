function [data] = fillStructLifetimeHist_segment(data, censor)
% fill experiment structure with lifetime histogram based on the lftInfo
% file saved at the specified directory; separate lftHist vectors are
% calculated based on the image segmentation in the specified folders
% 
% SYNOPSIS [exp]=fillStructLifetimeInfo(exp);
%
% INPUT     data: experiment structure, which has to contain the fields
%                .source
%                .movieLength
%                .segmentDataFileName 
%                .segmentDataFilePath 
%           censor: censor variable; if this is set to zero, all
%               trajectories regardless of status (complete or partial) are
%               counted for the lifetime histogram - default is censor=1
% OUTPUT    creates new fields in the data structure called
%                .lftHist   
%                .lftHist_InRegion
%                .lftHist_OutRegion
%           which contain the lifetime histogram of the objects 
%                 - in the entire image
%                 - only inside the segmented region
%                 - only outside the segmented region
% REMARKS 
%
% Dinah Loerke, last modified June 25, 2008
% Francois Aguet, 01/22/2010

if (nargin < 2)
    censor = 1;
end

for i = 1:length(data)
    
    fprintf('Movie no. %d\n', i);    
    load([data(i).source 'LifetimeInfo' filesep 'lftInfo.mat']);
    
    % check if segmentation status already exists
    if isfield(data(i), 'segmentStatus') && ~isempty(data(i).segmentStatus)
        segmentStatusVector = data(i).segmentStatus;
    else
        SegmentMask = imread([data(i).segmentDataFilePath data(i).segmentDataFileName]);
        % calculate segmentation status (1=Inside, 0=oustide segmented region)
        [segmentStatusVector] = calcIORegionLfthistSimple(lftInfo, SegmentMask);
        data(i).segmentStatus = segmentStatusVector';
    end
    
    lftMat = lftInfo.Mat_lifetime;
    statMat =  lftInfo.Mat_status;

    sx = size(lftMat, 1);

    lftVec = NaN(sx,1);
    lftVecIN = NaN(sx,1);
    lftVecOUT = NaN(sx,1);

    % DEFAULT:
    % a trajectory is counted for the lifetime analysis if the status of 
    % the trajectory is ==1 and the value of the gaps is ==4
    % IF censor==0
    % count all status trajectories
    for k=1:sx
        % current status vector
        cstat = nonzeros(statMat(k,:));
        % current lifetime vector
        clft = lftMat(k,:);
        
        countStat = ( (min(cstat)==1) & (max(cstat)<5) );
        if censor==0
            countStat = (max(cstat)<5);
        end
        
        if countStat==1
            lftVec(k) = max(clft);
            % count into either IN or OUT vector depending on segmentation status
            if segmentStatusVector(k)==1
                lftVecIN(k) = max(clft);
            else
                lftVecOUT(k) = max(clft);
            end
        end
    end   
    
    % lifetime vectors containing all counted lifetime lengths
    lftVec = lftVec(isfinite(lftVec));
    lftVecIN = lftVecIN(isfinite(lftVecIN));
    lftVecOUT = lftVecOUT(isfinite(lftVecOUT));

    % lifetime histograms
    data(i).lftHist = hist(lftVec, 1:data(i).movieLength);
    data(i).lftHist_InRegion = hist(lftVecIN, 1:data(i).movieLength);
    data(i).lftHist_OutRegion = hist(lftVecOUT, 1:data(i).movieLength);
end