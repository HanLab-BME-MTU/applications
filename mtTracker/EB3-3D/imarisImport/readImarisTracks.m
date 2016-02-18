function tracks=readImarisTracks(positionTextFile)
trackCell=csv2cell(positionTextFile);
%% Formating results
XYZ=cell2mat(trackCell(4:end,1:3));
TrackId=cell2mat(trackCell(4:end,8));
TrackId=TrackId-TrackId(1)+1;
Timing=cell2mat(trackCell(4:end,7));
featureIdx=zeros(size(Timing));
for t=1:max(Timing)
    fIdx=cumsum(Timing==t);
    featureIdx(Timing==t)=fIdx((Timing==t));
end
%%
tracks(max(TrackId))=TracksHandle();

for tIdx=1:max(TrackId)
    tIdxMask=(TrackId==tIdx);
    tracks(tIdx).x=XYZ(tIdxMask,1)'+1;
    tracks(tIdx).y=XYZ(tIdxMask,2)'+1;
    tracks(tIdx).z=XYZ(tIdxMask,3)'+1;
    tracks(tIdx).A=ones(size(XYZ(tIdxMask,3)))';
    tracks(tIdx).dx=0.5*ones(size(XYZ(tIdxMask,3)))';
    tracks(tIdx).dy=0.5*ones(size(XYZ(tIdxMask,3)))';
    tracks(tIdx).dz=0.5*ones(size(XYZ(tIdxMask,3)))';
    tracks(tIdx).dA=0.5*ones(size(XYZ(tIdxMask,3)))';
    tracks(tIdx).startFrame=min(Timing(tIdxMask));
    tracks(tIdx).endFrame=max(Timing(tIdxMask));
    tracks(tIdx).segmentStartFrame=tracks(tIdx).startFrame;
    tracks(tIdx).segmentEndFrame=tracks(tIdx).endFrame;
    tracks(tIdx).tracksFeatIndxCG=featureIdx(tIdxMask)';
end
