function idlist=calcXCoord(idlist,fieldname,centroid,rotMat)
%function to calculate the extremal coordinates from idlist for displaying

%SYNOPSIS idlist=calcXCoord(idlist,centroid,rotMat)
%
%INPUT idlist
%      fieldname: name of the new field in idlist(1).stats
%      centroid: list of centers by which coordinates are shifted
%      rotMat: matrix by which the coordinates are rotated (optional)
%
%OUTPUT idlist with new fields:
%           idlist(1).stats.fieldname
%           
%
%c: 1/03 Jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nFrames=size(idlist,2);
coord=[];


sCen=size(centroid);
switch (sCen(1)==nFrames)+(sCen(2)==nFrames)*2
    case 2
        centroid=centroid';
    case 0
        error('bad centroid list')
end
            

%check if rotMat
if nargin==3
    rotMat=ones(3,3,nFrames);
elseif nargin~=4
    error('wrong number of input arguments')
end

k=1;
%read coordinates
for t=1:nFrames
    if ~isempty(idlist(t).linklist)
        for i=1:max(idlist(t).linklist(:,2))
            coord=[coord;(rotMat(:,:,t)*(idlist(t).spot(i).coord-centroid(k,:))')'];
        end
        k=k+1;
    end
end

maxCoord=max(coord,[],1);
minCoord=min(coord,[],1);

eval(['idlist(1).stats.',fieldname,'=[minCoord;maxCoord];']);