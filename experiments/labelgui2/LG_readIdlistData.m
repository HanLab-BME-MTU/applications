function idlistData = LG_readIdlistData(idlist)
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
% there are 6 different flags at the moment: 3 spotFlags, 3 tagFlags
% add one more flag for deleted frame
flagList = zeros(nTimepoints,6);

for t = 1:nTimepoints
    if ~isempty(idlist(t).linklist)
        
        % remember good
        goodIdx(t) = 1;
        
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

% Remove superfluous columns
zeroCols = all(flagList==0,1);
if all(zeroCols)
    % keep at least 1 column
    zeroCols(1) = false;
end
flagList(:,zeroCols) = [];

idlistData.nSpots = nSpots;
idlistData.goodTimes = find(goodIdx);
idlistData.maxSpots = max(nSpots);
idlistData.maxTags = idlist(1).stats.maxColor;
idlistData.labelcolor = idlist(1).stats.labelcolor;
idlistData.labellist = idlist(1).stats.labellist;
idlistData.flagList = flagList;