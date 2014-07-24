function idlistData = LG_readIdlistData(idlist,dataProperties)
%LG_readIdlistData reads basic statistics from the idlist for labelgui2
% IdlistData has the fields:
% nSpots;
% goodTimes
% maxSpots
% maxTags
% labelcolor
% labellist
% flagList
% if fields are added here, update LG_navi_menuLoadIdlist

nTimepoints = length(idlist);
[nSpots,goodIdx] = deal(zeros(nTimepoints,1));
% there are 6 different flags at the moment: 3 spotFlags, 3 tagFlags,
% one for close to border
% add one more flag for deleted frame
flagList = zeros(nTimepoints,6);

for t = 1:nTimepoints
    if ~isempty(idlist(t).linklist)
        
        % remember good
        goodIdx(t) = true;
        
        % count nSpots. Account for fusions
        nSpots(t) = ...
            nnz(idlist(t).linklist(:,2) > 0 & idlist(t).linklist(:,3)~=4);
        
        % collect flags. 
        flags1 = idlist(t).linklist(:,3);
        flags2 = idlist(t).linklist(:,5);
        % flags 3 are very much not interesting - don't store
        f3 = flags2 == 3;
        flags1(f3) = [];
        flags2(f3 | flags2==0) = [];
        
        flags = unique([flags1;10+flags2]);
        flagList(t,1:length(flags)) = flags';        
    else
        flagList(t,1) = 21;
    end
end

goodTimes = find(goodIdx);

% flag 22: position close to border (within two psf-widths)
lowerBorder = 1.5*dataProperties.FT_SIGMA.*[dataProperties.PIXELSIZE_XY,...
    dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_Z];
upperBorder = [dataProperties.PIXELSIZE_XY, dataProperties.PIXELSIZE_XY,...
    dataProperties.PIXELSIZE_Z].* (dataProperties.movieSize(1:3)+1)...
    - lowerBorder;
pos = catStruct(3,'idlist.linklist(:,9:11)');
closeIdx = any(pos(:,1,:) < lowerBorder(1) | pos(:,1,:) > upperBorder(1) | ...
    pos(:,2,:) < lowerBorder(2) | pos(:,2,:) > upperBorder(2) | ...
    pos(:,3,:) < lowerBorder(3) | pos(:,3,:) > upperBorder(3),2);
% check for single occurences - estimated single occurences shouldn't count
% as close to border!
goodCheck = catStruct(3,'idlist.linklist(:,5)');
badIdx = goodCheck == 3;
closeIdx(badIdx) = false;
closeIdx = any(closeIdx,1);
flagList(goodTimes(closeIdx),end+1) = 22;


% Remove superfluous columns
zeroCols = all(flagList==0,1);
if all(zeroCols)
    % keep at least 1 column
    zeroCols(1) = false;
end
flagList(:,zeroCols) = [];

idlistData.nSpots = nSpots;
idlistData.goodTimes = goodTimes;
idlistData.maxSpots = max(nSpots);
idlistData.maxTags = idlist(1).stats.maxColor;
idlistData.labelcolor = idlist(1).stats.labelcolor;
idlistData.labellist = idlist(1).stats.labellist;
idlistData.flagList = flagList;