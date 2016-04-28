function tracks=readImarisTracks(positionTextFile,varargin)
ip=inputParser();
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('positionTextFile')
ip.addOptional('intensityFile',[],@(x) ischar(x));
ip.parse(positionTextFile,varargin{:})

%trackCell=csv2cell(positionTextFile);
posTab = readtable(positionTextFile,'Format','%f%f%f%s%s%s%f%f%f%[^\n\r]','HeaderLines',3);
trackCell=table2cell(posTab);
if(~isempty(ip.Results.intensityFile))
    intensityTab = readtable(ip.Results.intensityFile,'Format','%f%s%s%f%f%f%f%[^\n\r]','HeaderLines',3);
    intCell=table2cell(intensityTab);
end

%% Formating results
TrackIdCell=trackCell(:,8);
correctLine=cellfun(@(x) isnumeric(x),TrackIdCell);

XYZ=cell2mat(trackCell(correctLine,1:3));
A=cell2mat(intCell(correctLine,1));
TrackId=cell2mat(trackCell(correctLine,8));
TrackId=TrackId-TrackId(1)+1;
Timing=cell2mat(trackCell(correctLine,7));
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
    tracks(tIdx).A=A(tIdxMask)';
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
