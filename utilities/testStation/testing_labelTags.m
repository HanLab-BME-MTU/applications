function idlist = testing_labelTags(idlist,positions,dataProperties)
%TESTING_LABELTAGS labels tags in idlists according to the ground truth
%
% SYNOPSIS: idlist = testing_labelTags(idlist,positions)
%
% INPUT idlist: idlist
%		positions: ground truth. nTimepoints-by-4-by-nTags array
%		(amplitudes aren't being used) 
%       dataProperties
%
% OUTPUT idlist: idlist with labelcolor according to nTag.
%
% REMARKS - Only first timepoint is being used!
%         - The code swaps x,y of positions to make the correspondence
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: jdorn
% DATE: 27-May-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find first good frame and read positions. Read everything, including
% estimates!
linklists = cat(3,idlist.linklist);
linkerNTags = size(linklists,1);
linkerPositions = linklists(:,9:11,1)./...
    repmat([dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_Z],linkerNTags,1);

goodTimes = squeeze(linklists(1,1,:));


% read positions from position-structure (swap xy) make nTags-by-3 array
groundTruth = positions(goodTimes(1),[2,1,3],:);
groundTruth = permute(groundTruth,[3,2,1]);
nTags = size(positions,3);

% LAP. Get index into ground truth for every linked tag. Set maximum
% allowed deviation to three pixels.  
distanceMatrix = distMat2(linkerPositions,groundTruth);
[tagNumber] = lap(distanceMatrix,-1,0,1,3); 

% loop through tags. Add numbers, set good/bad
for i = 1:linkerNTags
    if tagNumber(i) > nTags
        % bad tag
        idlist = LG_setGoodTag(idlist,goodTimes,i,dataProperties,0);
    else
        % good tag
        idlist = LG_setGoodTag(idlist,goodTimes,i,dataProperties,1);
        idlist(1).stats.labelcolor{i} = num2str(tagNumber(i));
    end
end
