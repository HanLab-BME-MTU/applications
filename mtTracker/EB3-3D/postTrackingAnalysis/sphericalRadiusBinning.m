function [sphCumulCoordCell]=sphericalRadiusBinning(sphCoord,radius,timeInterval,tracks,poleId)
%% Bin the spherical coordinate function of their radius. Optionnaly provides tracking ID too. 
%  sphCoord is a structure defining .azimuth{fIdx}, .elevation{fIdx} and
%  .rho{fIdx}
%  radius is a matrix defining the binning edges. 
% 

% Concatenate the per-frame data-structure
cumulAzi=vertcat(sphCoord.azimuth{:}); 
cumulElev=vertcat(sphCoord.elevation{:});
cumulRho=vertcat(sphCoord.rho{:});
if(~isempty(poleId))
    cumulPoleId=vertcat(poleId{:});
end
T=arrayfun(@(f) ones(length(sphCoord.azimuth{f}),1)*f*timeInterval,(1:length(sphCoord.azimuth)),'unif',0);
cumulTimePt=vertcat(T{:});

% Associate each detection with their a trackID
trackIdCell=cellfun(@(x) ones(size(x)),sphCoord.azimuth,'unif',0);
for tIdx=1:length(tracks)
    aTrack=tracks(tIdx);
    notGapIdx=find(~aTrack.gapMask);
    F=aTrack.startFrame:aTrack.endFrame;
    for tpIdx=notGapIdx
        trackIdCell{F(tpIdx)}(aTrack.tracksFeatIndxCG(tpIdx))=tIdx;
    end   
end
cumulTrackId=vertcat(trackIdCell{:});

% Bin each Kinetochore according to <radius>
minKinRho=min(cumulRho,[],2);
[~,binIdx]=histc(minKinRho,radius);
sphCumulCoordCell=cell(1,length(radius));   
for rIdx=1:length(radius)
    sphCumulCoordCell{rIdx}.azimuth=cumulAzi(binIdx==rIdx,:);
    sphCumulCoordCell{rIdx}.elevation=cumulElev(binIdx==rIdx,:);
    sphCumulCoordCell{rIdx}.rho=cumulRho(binIdx==rIdx,:);
    if(~isempty(poleId))
        sphCumulCoordCell{rIdx}.poleId=cumulPoleId(binIdx==rIdx,:);
    end
    sphCumulCoordCell{rIdx}.time=cumulTimePt(binIdx==rIdx);
    sphCumulCoordCell{rIdx}.trackId=cumulTrackId(binIdx==rIdx);
end
