% Augment the structures with spherical Coordinate. 
%% Load the pole info
outputDirPoleDetect='Y:\bioinformatics\Danuser_lab\externBetzig\analysis\proudot\anaProject\sphericalProjection\prometaphase\cell1_12_half_volume_double_time\EB3\poles\simplex_scale_003'
tmp=load([outputDirPoleDetect filesep 'poleDetection.mat']);
poleMovieInfo=tmp.poleMovieInfo;

% WARNING: this is not a trajectory, merely a collection of poles to ease
% implementation.
P1=struct();
P1.x=arrayfun(@(d) MD.pixelSize_*(d.xCoord(1,1)-1)+1,poleMovieInfo)';
P1.y=arrayfun(@(d) MD.pixelSize_*(d.yCoord(1,1)-1)+1,poleMovieInfo)';
P1.z=arrayfun(@(d) MD.pixelSize_*(d.zCoord(1,1)-1)+1,poleMovieInfo)';

P2=struct();
P2.x=arrayfun(@(d) MD.pixelSize_*(d.xCoord(2,1)-1)+1,poleMovieInfo)';
P2.y=arrayfun(@(d) MD.pixelSize_*(d.yCoord(2,1)-1)+1,poleMovieInfo)';
P2.z=arrayfun(@(d) MD.pixelSize_*(d.zCoord(2,1)-1)+1,poleMovieInfo)';

refP1=FrameOfRef();
refP1.setOriginFromTrack(P1);
refP1.setZFromTrack(P2);
refP1.genBaseFromZ();

refP2=FrameOfRef();
refP2.setOriginFromTrack(P2);
refP2.setZFromTrack(P1);
refP2.genBaseFromZ();

%%
tmp=load('Y:\bioinformatics\Danuser_lab\externBetzig\analysis\proudot\anaProject\sphericalProjection\prometaphase\cell1_12_half_volume_double_time\Kin\track\tracksLabRef.mat')
tracks=tmp.tracksLabRef
couilleIdx=zeros(1,length(tracks));
xmag=zeros(1,length(tracks));
ymag=zeros(1,length(tracks));

for i=1:length(tracks)
tr=tracks(i);
refP1K=FrameOfRef;
refP1K.origin=refP1.origin(tr.f,:);
refP1K.setZFromTrack(tr);
refP1K.genBaseFromZ;
refP1K.applyBaseToTrack(tr,'P1K');
    if(max(abs(tr.P1K.x))>1)
        couilleIdx(i)=1;
    end
    xmag(i)=max(abs(tr.P1K.x));
    ymag(i)=max(abs(tr.P1K.y));

end
