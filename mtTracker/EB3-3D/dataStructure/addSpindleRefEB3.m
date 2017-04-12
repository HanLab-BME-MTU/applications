function [EB3Tracks]=addSpindleRefEB3(MD,poleMovieInfo,EB3Tracks,EB3SphCoord,EB3Inliers,EB3PoleId)

 % WARNING: this is not a trajectory, merely a collection of poles to ease
% implementation.
P1=struct();
P1.x=arrayfun(@(d) MD.pixelSize_*(d.xCoord(1,1)-1)+1,poleMovieInfo)';
P1.y=arrayfun(@(d) MD.pixelSize_*(d.yCoord(1,1)-1)+1,poleMovieInfo)';
P1.z=arrayfun(@(d) MD.pixelSize_*(d.zCoord(1,1)-1)+1,poleMovieInfo)';
P1.f=1:length(poleMovieInfo);

P2=struct();
P2.x=arrayfun(@(d) MD.pixelSize_*(d.xCoord(2,1)-1)+1,poleMovieInfo)';
P2.y=arrayfun(@(d) MD.pixelSize_*(d.yCoord(2,1)-1)+1,poleMovieInfo)';
P2.z=arrayfun(@(d) MD.pixelSize_*(d.zCoord(2,1)-1)+1,poleMovieInfo)';
P2.f=1:length(poleMovieInfo);


refP1=FrameOfRef();
refP1.setOriginFromTrack(P1);
refP1.setZFromTrack(P2);
refP1.genBaseFromZ();

refP2=FrameOfRef();
refP2.setOriginFromTrack(P2);
refP2.setZFromTrack(P1);
refP2.genBaseFromZ();

poleRefs=[refP1 refP2];

% For MT detection
% tic;
% detLabRef=Detections(EB3LabRef);
% detLabRef.scale(MD);
% dp1=refP1.applyBaseToDetection(detLabRef,'pole1');
% dp2=refP2.applyBaseToDetection(detLabRef,'pole2');
% dp1.addSphericalCoord();
% dp2.addSphericalCoord();
% toc

%% For MT tracks
tic;
for tIdx=1:length(EB3Tracks)
    %progressText(tIdx/length(EB3tracks),'Loading EB3 spherical coordinates.');

    tr=EB3Tracks(tIdx);
     try
        tr.addprop('inliers');
        tr.addprop('poleId');
        tr.addprop('azimuth');      % DEPRECATED
        tr.addprop('elevation');    % DEPRECATED
        tr.addprop('rho');          % DEPRECATED
    catch
    end;
    tr.x=(tr.x-1)*MD.pixelSize_+1;
    tr.y=(tr.y-1)*MD.pixelSize_+1;
    tr.z=(tr.z-1)*MD.pixelSize_+1;

    nonGap=~tr.gapMask();
    tr.poleId=nan(size(tr.f));
    tr.poleId(nonGap)=arrayfun(@(i,f) EB3PoleId{f}(i), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
    tr.inliers=nan(size(tr.f));
    tr.inliers(nonGap)=arrayfun(@(i,f) EB3Inliers{f}(i), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));

    % DEPRECATED
    tr.azimuth=nan(2,length(tr.f));
    tr.elevation=nan(2,length(tr.f));
    tr.rho=nan(2,length(tr.f));
    for poleID=1:2
        tr.azimuth(poleID,nonGap)=arrayfun(@(i,f) EB3SphCoord.azimuth{f}(i,poleID), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
        tr.elevation(poleID,nonGap)=arrayfun(@(i,f) EB3SphCoord.elevation{f}(i,poleID), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
        tr.rho(poleID,nonGap)=arrayfun(@(i,f) EB3SphCoord.rho{f}(i,poleID), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
    end
    % END DEPRECATED
end
toc;
%%
tic;
refP1.applyBaseToTrack(EB3Tracks,'pole1');
refP2.applyBaseToTrack(EB3Tracks,'pole2');
%%
% augment pole ref projected data with spherical coordinate
refName={'pole1','pole2'};
for tIdx=1:length(EB3Tracks)
%     trackPoleRefs=[];
    for poleID=1:length(poleRefs);
        % Copying EB3 track
        tr=getfield(EB3Tracks(tIdx),refName{poleID});
        % Adding correponding spherical coordinate
        try
            tr.addprop('azimuth');
            tr.addprop('elevation');
            tr.addprop('rho');
        catch
        end;

        nonGap=~tr.gapMask();
        tr.azimuth=nan(1,length(tr.f));
        tr.elevation=nan(1,length(tr.f));
        tr.rho=nan(1,length(tr.f));

        tr.azimuth(nonGap)=arrayfun(@(i,f) EB3SphCoord.azimuth{f}(i,poleID), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
        tr.elevation(nonGap)=arrayfun(@(i,f) EB3SphCoord.elevation{f}(i,poleID), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
        tr.rho(nonGap)=arrayfun(@(i,f) EB3SphCoord.rho{f}(i,poleID), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
    end
end
