function [data] = fillStructDensities_segment(data,censor,outvar)
% fill experiment structure with density values for the inside and outside
% regions based on the image segmentation in the specified folders
% 
% SYNOPSIS [data] = fillStructDensities_segment(data,censor,outvar);
%
% INPUT     data: experiment structure, which has to contain the fields
%                .source
%                .movieLength
%                .segmentDataFileName 
%                .segmentDataFilePath 
%                AND, if outvar==1, additionally the fields
%                .segmentDataFilePathOUT
%                .segmentDataFileNameOUT
%           censor: censor variable (optional)
%                   if censor==0, all trajectories regardless of status 
%                   (complete or partial) are counted for the object
%                   density
%                   DEFAULT is censor==1
%           outvar : outside segmentation variable (optional)
%                   if outvar==1, then the densities in the outside region
%                   are also restricted, to the segmented region specified
%                   in the fields .segmentDataFileNameOUT - this segmented
%                   region could e.g. contain the segmentation of the cell
%                   outline, as opposed to the inside segmentation of cell
%                   subregions
% OUTPUT    creates new fields in the data structure called
%                .density   
%                .density_InRegion
%                .density_OutRegion
%           which contain the object densities 
%                 - in the entire image
%                 - only inside the segmented region
%                 - only outside the segmented region
%           
%           NOTE1: the dimension of the density is 
%           # of objects per square pixel
%
%       NOTE2: Unlike in the function determinePitDensities, the densities 
%       here are calculated based NOT only on the number of all detected 
%       objects in every frame, but using the full trajectories (i.e. the 
%       trajectories that are also counted towards the survival function).
%       Trajectories can be excluded if they have the wrong 'censor'
%       status or if they aren't counted towards the lifetime analysis
%       (status=5). In addition, an object's segmentation status (i.e.
%       whether it's considered to be inside or outside the pattern) is 
%       based on its entire trajectory, and not only its current position.
%       This approach can slightly reduce the number of counted objects 
%       compared to the function determinePitDensities. 
%
% Dinah Loerke, modified July 23, 2008
% Dinah Loerke, modified July 29, 2008
% Francois Aguet, Feb 2010

if (nargin < 2)
    censor = 1;
end
if (nargin < 3)
    outvar = 0;
end

for i = 1:length(data)
    
    fprintf('Movie no. %d\n', i);
    load([data(i).source 'LifetimeInfo' filesep 'lftInfo.mat']);
    
    % check if segmentation status already exists
    if isfield(data,'segmentStatus')
        segmentStatusVector = data(i).segmentStatus;
    else
        % load segmentation image from specified location
        SegmentMask = imread([data(i).segmentDataFilePath data(i).segmentDataFileName]);    
        if ~islogical(SegmentMask)
            SegmentMask = logical(SegmentMask/max(SegmentMask(:)));
        end
        % calculate segmentation status (1=Inside, 0=oustide segmented region)
        [segmentStatusVector] = calcIORegionLfthistSimple(lftInfo, SegmentMask);
    end
    
    % if outvar==1, load segmentation data to restrict the OUTSIDE area
    % additionally (e.g. to the cell outline)
    if outvar==1
        % load segmentation image from specified location
        SegmentMaskOUT = imread([data(i).segmentDataFilePathOUT data(i).segmentDataFileNameOUT]);
        if ~islogical(SegmentMaskOUT)
            SegmentMaskOUT = logical(SegmentMaskOUT/max(SegmentMaskOUT(:)));
        end
    end
    
    % calculate average number of inside and outside objects for all
    % appropriate objects (depending on censor status)
    
    lftMat = lftInfo.Mat_lifetime;
    statMat =  lftInfo.Mat_status;

    sx = size(lftMat, 1);

    % == DEFAULT:
    % a trajectory is counted for the density analysis if the status of 
    % the trajectory is ==1 and the value of the gaps is ==4
    % == IF censor==0
    % count all status trajectories
    
    countMat = 0 * full(statMat);
    
    for k=1:sx
        % current status vector
        cstat = nonzeros(statMat(k,:));
                
        countStat = ( (min(cstat)==1) & (max(cstat)<5) );
        if censor==0
            countStat = (max(cstat)<5);
        end
        
        % if countStat==1, i.e. if this trajectory fits the censoring
        % requirements, then all the detected points in the trajectory are
        % counted
        if countStat==1
            countMat(k,:) = (statMat(k,:)<4) & (statMat(k,:)>0);
        end
    end
    
    % split countMat into In and OUT matrices depending on segmentation 
    % status determined by segmentStatusVector
    countMatIN = countMat(segmentStatusVector==1,:);
    countMatOUT = countMat(segmentStatusVector~=1,:);
            
    % number of points inside and outside, averaged over all frames in the movie
    numIN = nanmean(sum(countMatIN,1));
    numOUT = nanmean(sum(countMatOUT,1));
    numALL = nanmean(sum(countMat,1));    
        
    % calculate density by dividing by segmented area - depending on
    % outvar, restrict outside region to segmented image as well
        
    if outvar==1
        areaIN = sum(SegmentMaskOUT(:) & SegmentMask(:));
        areaOUT = sum(SegmentMaskOUT(:) & ~SegmentMask(:));
        areaALL = sum(SegmentMaskOUT(:));
    else
        areaIN = sum(SegmentMask(:));
        areaOUT = sum(~SegmentMask(:));
        areaALL = length(SegmentMask(:));
    end
    
    densityIN = numIN/areaIN;
    densityOUT = numOUT/areaOUT;
    densityALL = numALL/areaALL;

    % lifetime histograms
    data(i).density = densityALL;
    data(i).density_InRegion = densityIN;
    data(i).density_OutRegion = densityOUT;
end