function region=maskSpots(sl,fSze,dataProperties)
%MASKSPOTS craete a region struct containing mask for a single frame
%
% SYNOPSIS  region=maskSpots(sl,fSze)
%
% INPUT     sl:    list of spots from spotfinder (2 frames)
%           fSze:  full size of frame
% 
% OUTPUT nsl : tracked frame spots info
%   
% c: 07/08/02	dT

% Constants

%global FT_SIGMA PIXELSIZE_XY PIXELSIZE_Z MAXSPOTS;
if nargin <3 | isempty(dataProperties)
    %load_sp_constants;
    error('Missing dataProperties');
else
    
    PIXELSIZE_XY = dataProperties.PIXELSIZE_XY;
    PIXELSIZE_Z = dataProperties.PIXELSIZE_Z;
    FT_SIGMA = dataProperties.FT_SIGMA;
    MAXSPOTS = dataProperties.MAXSPOTS;
end

%pixel <-> micron
p2m=ones(MAXSPOTS,1)*[PIXELSIZE_XY PIXELSIZE_XY PIXELSIZE_Z];


%number of spots from spfinder
nspCur=max(sl(1).linklist(:,2));

[tval,uniqSpotIdx,invUniqSpotIdx] = unique(sl(1).linklist(:,2));

%coordinate in pixels (matlab convention x<->y)
coordsCurMu=sl(1).linklist(uniqSpotIdx,9:11);
coordsCurPix=[coordsCurMu(:,2) coordsCurMu(:,1) coordsCurMu(:,3)]./p2m(1:nspCur,:);

%reconstruct spot amplitude
for i=1:max(sl(1).linklist(:,2))
    rowIdx=find(sl(1).linklist(:,2)==i);
    spotAmp(i)=sum(sl(1).linklist(rowIdx,8));
end

maxDynRange=max(spotAmp);

region=[];
% mask spots
for s=1:nspCur
    region(s).center=coordsCurPix(s,:);
    region(s).amp=spotAmp(sl(1).linklist(uniqSpotIdx(s),2));
    gm=gaussMaskThres(region(s).amp,FT_SIGMA,fSze,coordsCurPix(s,:),maxDynRange/50);
    region(s).coords=gm(:,1);
    region(s).gauss=gm(:,2);
    %check for overlap with previous masks
    region(s).ovlp=[s];
    for i=1:length(region)-1
        if ~isempty(intersect(region(i).coords,region(s).coords))
            region(i).ovlp=[region(i).ovlp s];
            region(s).ovlp=[region(s).ovlp i];
        end;
    end;
end;

%  calc for overlaps
for s=1:nspCur
    tempVal=region(s).gauss;
    for ov=region(s).ovlp(2:end)
        [its ia ib]=intersect(region(s).coords, region(ov).coords);
        tempVal(ia)=tempVal(ia)+region(ov).gauss(ib);
    end;
    region(s).ratio=region(s).gauss./tempVal;
end;

%back to 3d coords
for s=1:nspCur
[x y z]=ind2sub(fSze,region(s).coords);
region(s).coords=[x y z];
end;