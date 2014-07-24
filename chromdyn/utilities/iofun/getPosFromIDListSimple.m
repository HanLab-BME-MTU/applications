function tagCoord = getPosFromIDListSimple(idlist)

numTP = length(idlist);

numTags = size(idlist(1).linklist,1);

tagCoord = NaN(numTags,3,numTP);

for iTP = 1 : numTP
    if ~isnan(idlist(iTP).linklist)
        tagCoord(:,:,iTP) = idlist(iTP).linklist(:,9:11);
    end
end
