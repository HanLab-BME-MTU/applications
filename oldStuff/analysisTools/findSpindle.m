function [spbList,cenList]=findSpindle(idlist)
%returns tag colors of spb/cen in idlist
%
%SYNOPSIS [spbList,cenList]=findSpindle(idlist)
%
%INPUT idlist with labeled tags
%
%OUTPUT list of tag colors of spb (2) and cen (0-2)
%
%c: 1/03 Jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nColors=size(idlist(1).linklist,1);
spbList=[];
cenList=[];
for i=1:nColors
    if findstr(idlist(1).stats.labelcolor{i},idlist(1).stats.labellist{1})
        spbList=[spbList;i];
    else
        cenList=[cenList;i];
    end
end