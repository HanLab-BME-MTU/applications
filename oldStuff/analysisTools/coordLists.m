function varargout=coordLists(idlist, tagList, nFrames)
%coordLists extracts spb/cen coordinates from idlist
%
%SYNOPSIS varargout=coordLists(idlist, tagList, nFrames)
%
%INPUT idlist 
%      tagList: list of tag colors created with findSpindle
%
%OUTPUT nTags nFramesx3-matrices containing the coordinates labeled with
%       the respective tags
%
%c: 1/03 Jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
out=zeros(nFrames,3);

for i=1:nargout
    for j=1:nFrames
        if ~isempty(idlist(j).linklist)
            rowidx=find(idlist(j).linklist(:,4)==2^(tagList(i)-1));
            spotnum=idlist(j).linklist(rowidx,2);
            coord=idlist(j).spot(spotnum).coord;
            out(j,:)=coord;
        end
    end
    varargout(i) = {out};
end