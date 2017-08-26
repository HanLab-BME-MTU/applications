function [EB3Tracks]=addSpindleRefEB3(MD,poleRefs,EB3Tracks,EB3SphCoord,EB3Inliers,EB3PoleId)

% For MT detection
% tic;
% detLabRef=Detections(EB3LabRef);
% detLabRef.scale(MD);
% dp1=refP1.applyBaseToDetection(detLabRef,'pole1');
% dp2=refP2.applyBaseToDetection(detLabRef,'pole2');
% dp1.addSphericalCoord();
% dp2.addSphericalCoord();
% toc
refP1=poleRefs(1);
refP2=poleRefs(2);

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
