function [pixelsize, extendedStats] = calculatePixelSize(image,linetype,trueDistance,trueDistanceType)
%CALCULATEPIXELSIZE calculates the pixel size from a regular grid of squares
%
% SYNOPSIS [pixelsize, extendedStats] = calculatePixelSize(image,linetype,trueDistance,trueDistanceType)
%
% INPUT    image:            grayvalue image of regular grid of empty squares
%          trueDistance:     known distance
%          trueDistanceType: 1: distance is measured from left (right) edge 
%                               of square to left (right) edge of adjacent square 
%          linetype:         [1/{2}] 1: bright line on blackground, 2: black
%                               line on white background
%
% OUTPUT   pixelsize:        um/pixel
%          extendedStats:    cell with individual distances, angle
%
% The grid should be roughly aligned with the image.
% If the lineSpacing graphs show include too much, put tolerance up
%
%c: 04/04 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% no test input - we assume the user knows what he/she's doing

% tolerance defines down to how many percent of the maximum response of the
% hough transform is used to calculate the center of a line
tolerance = 0.1;

%===============
% EXTRACT LINES
%===============

% prepare image
imgSize = size(image);

% convert to grayscale if necessary
if length(imgSize) == 3 & imgSize(3) == 3
    image = sum(image,3);
end

% scale to [0,1]
image = image - min(image(:));
image = image/max(image(:));

% prepare options for imLineExtract
opt.scales = [1:5:101];
opt.linetype = linetype;

% extract lines
[binMask,ori] = imLineExtract(image,opt);

% cleanup maps (not done anymore)

% bwlabel
% labeledResp = bwlabel(binMask);
% [labels,numOccur,label2map] = countEntries(labeledResp);
% numOccur(1) = 0; % we do not care about all that's not a line
% 
% % use robust mean to find both the group with few and with a lot of pixels
% [mean1,std1,include1] = robustMean(numOccur);
% exclude1 = missingIndices(include1,length(numOccur));
% 
% if ~isempty(exclude1) % maybe there are no small groups
%     [mean2,std2,include2] = robustMean(numOccur(exclude1));
%     include2 = exclude1(include2);
%     exclude2 = missingIndices(include2,length(numOccur));
% else
%     mean2 = -999;
%     exclude2 = [1:length(numOccur)];
% end
% 
% % normally, the group with the larger mean should be the group that we are
% % interested in.
% % goodGroup are the indices into resp, ori that correspond to large groups
% % careful! labeledResp starts at 0!
% if mean1>mean2
%     goodGroupIdx = find(ismember(labeledResp(:)+1,include1));
%     badGroupIdx  = find(ismember(labeledResp(:)+1,exclude1));
% elseif mean1<mean2
%     goodGroupIdx = find(ismember(labeledResp(:)+1,include2));
%     badGroupIdx  = find(ismember(labeledResp(:)+1,exclude2));
% else
%     error('can''t separate large groups from small ones')
% end
% 
% % cleanup resp
% binMask(badGroupIdx) = 0;
% 
% % cleanup ori. Angle goes from 0 to 2pi, everything else will be -1
% ori(find(binMask==0)) = -1;

% show cleaned binMask and image
uivpH = uiviewpanel;
set(uivpH,'Name','raw image')
imshow(image,[]);

uivpH = uiviewpanel;
set(uivpH,'Name','extracted lines')
imshow(binMask);


%==================


%=====================
% GET LINE INFO
%=====================

% we now perform a Hough transform on the image. There should be peaks
% around 0 degree and 90 degree; we find the angles and then go and measure
% distances

% Hough transform
theta = [-45:134];
[hough,xp] = radon(binMask,theta);

% show result
figure('Name','HoughTransform'), imagesc(theta, xp, hough); colormap(hot);

% we have two kinds of lines, so we want to run this twice
for iLin = 1:2
% the maxima should be along the same angle
[maxResp,maxIdx] = max(hough(:));
[posIdx,angIdx] = ind2sub(size(hough),maxIdx);

% now find the distances associated with the angle
distMap = hough(:,angIdx);

% find clusters along the line

% make 0/1-map, taking only entries with more than 5% of lineMax
    zeroOneList = zeros(size(distMap));
    zeroOneList(find(distMap>tolerance*maxResp)) = 1;
    
    figure('Name',['Line Spacing for angle ' num2str(theta(angIdx))]),plot(distMap),hold on,plot(zeroOneList*maxResp,'-r')
    
    distGroupIdx = findGroups(zeroOneList(:),1);
    
    % join groups that are separated by almost nothing
    groupDist = distGroupIdx(2:end,1) - distGroupIdx(1:end-1,2);
    joinIdx = find(groupDist < 4);
    if ~isempty(joinIdx)
        % make end of second group end of first group
        distGroupIdx(joinIdx,2) = distGroupIdx(joinIdx+1,2);
        % remove other group
        distGroupIdx(joinIdx+1,:) = [];
    end
    
    % now we can loop throught the groups and find their weighted
    % center-index
    for iGroup = 1:size(distGroupIdx,1);
        dmIdx = [distGroupIdx(iGroup,1):distGroupIdx(iGroup,2)]';
        weightedCenter(iGroup) = weightedStats(dmIdx,distMap(dmIdx),'w');
    end
    
    % store data
    lineCenterData{iLin,1} = theta(angIdx);
    lineCenterData{iLin,2} = weightedCenter;
    
    % mask houghTransform
    hough(:,angIdx-10:angIdx+10) = 0;
end

%===========================
% CALCULATE PIXELSIZE
%===========================

switch trueDistanceType
    case 1 % API - standard
        
        % the distance is measured from the left (right) edge of a square to the
        % left (right) edge of the adjacent square
        for i = 1:2
            wc = lineCenterData{i,2};
            deltas = wc(3:end) - wc(1:end-2);
            
            % we average deltas, then divide (MUCH better)
            avgPixelsize(i,1) = trueDistance/mean(deltas);
            avgPixelsize(i,2) = length(deltas);
            
        end
        
        % now we average over the two directions
        pixelsize = sum(prod(avgPixelsize,2))/sum(avgPixelsize(:,2));
        
        extendedStats = lineCenterData;
        
    otherwise
        error('this option is not implemented yet')
end