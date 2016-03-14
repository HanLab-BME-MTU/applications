function [azimuth,elevation,time,trackId,poleId]=sphericalDistribution(tracks,azimuthCells,elevationCells,rhoCells,poleIdCells,radius)
% Todo:  preeallocate and strip. 

crossingEB3 = repmat(struct('xCoord', [], 'yCoord',[],'zCoord', [], 'amp', [], 'int',[]),max([tracks.endFrame]),1);
interpCrossingEB3=repmat(struct('xCoord', [], 'yCoord',[],'zCoord', [], 'amp', [], 'int',[]),max([tracks.endFrame]),1);
cumulInterpCrossingEB3=interpCrossingEB3(1);
azimuth=[];
elevation=[];
time=[];
trackId=[];
poleId=[];


%% allocation of Amira coloring properties 
% detections
inOutProperties=cell(max([tracks.endFrame]),1);
association=cell(max([tracks.endFrame]),1);
% tracks
crossingTracks=zeros(length(tracks),1);

% GAP filling using the last known position
se=[zeros(1,tracks.numTimePoints) 1 ones(1,tracks.numTimePoints)];

for tIdx=1:length(tracks)
    aTrack=tracks(tIdx);
    gi=aTrack.gapMask;
    F=aTrack.f;
    if(any(gi))
        copyIdx=1:aTrack.lifetime;
        copyIdx(gi)=0;
        copyIdx=imdilate(copyIdx,se);
        aTrack.x=aTrack.x(copyIdx);
        aTrack.y=aTrack.y(copyIdx);
        aTrack.z=aTrack.z(copyIdx);
        aTrack.tracksFeatIndxCG=aTrack.tracksFeatIndxCG(copyIdx);
        F=aTrack.f(copyIdx);
    end
    

    %Collect rho,azimuth and elevation values
    tracksRhos=arrayfun(@(f,i) rhoCells{f}(i) ,F,aTrack.tracksFeatIndxCG);
    tracksAzimuth=arrayfun(@(f,i) azimuthCells{f}(i),F,aTrack.tracksFeatIndxCG);
    tracksElevation=arrayfun(@(f,i) elevationCells{f}(i),F,aTrack.tracksFeatIndxCG);
    tracksPoleId=arrayfun(@(f,i) poleIdCells{f}(i),F,aTrack.tracksFeatIndxCG);
    % find the first detection over the radius
    crossIdx=find(tracksRhos>radius,1);
    % If there is detection over the radius and 
    % if the first detection over the raidus is not the first detection of the track 
    if(~isempty(crossIdx))&&(crossIdx>1) 
        % timing
        frameCrossing=F(crossIdx-1);
        frameCrossed=F(crossIdx);
        
        % Interpolate Amizuth and Elevation of crossing and convert back to
        % euclidian referential
        rhoCrossing=tracksRhos(crossIdx-1); rhoCrossed=tracksRhos(crossIdx);
        azimuthCrossing=tracksAzimuth(crossIdx-1); azimuthCrossed=tracksAzimuth(crossIdx);
        elevationCrossing=tracksElevation(crossIdx-1); elevationCrossed=tracksElevation(crossIdx);
                        
        rhoRatio=(radius-rhoCrossing)/(rhoCrossed-rhoCrossing);
        azimuth=[azimuth azimuthCrossing+rhoRatio*(azimuthCrossed-azimuthCrossing)];
        elevation=[elevation elevationCrossing+rhoRatio*(elevationCrossed-elevationCrossing)];
        time=[time aTrack.t(crossIdx-1)+rhoRatio*(aTrack.t(crossIdx)-aTrack.t(crossIdx-1))];
        trackId=[trackId tIdx];
        poleId=[poleId tracksPoleId(1)];
        
        % Cartesian coord. for amira    
        xCross=aTrack.x(crossIdx-1)+rhoRatio*(aTrack.x(crossIdx)-aTrack.x(crossIdx-1));
        yCross=aTrack.y(crossIdx-1)+rhoRatio*(aTrack.y(crossIdx)-aTrack.y(crossIdx-1));
        zCross=aTrack.z(crossIdx-1)+rhoRatio*(aTrack.z(crossIdx)-aTrack.z(crossIdx-1));
        for j=[frameCrossing frameCrossed]
            interpCrossingEB3(j).xCoord=[interpCrossingEB3(j).xCoord; [xCross 0.5]];
            interpCrossingEB3(j).yCoord=[interpCrossingEB3(j).yCoord; [yCross 0.5]];
            interpCrossingEB3(j).zCoord=[interpCrossingEB3(j).zCoord; [zCross 0.5]];
        end
        
        % Filling structures for Amira
        crossingTracks(tIdx)=1;
        crossingEB3(frameCrossing).xCoord=[crossingEB3(frameCrossing).xCoord; [aTrack.x(crossIdx-1) 0.5]];
        crossingEB3(frameCrossing).yCoord=[crossingEB3(frameCrossing).yCoord; [aTrack.y(crossIdx-1) 0.5]];
        crossingEB3(frameCrossing).zCoord=[crossingEB3(frameCrossing).zCoord; [aTrack.z(crossIdx-1) 0.5]];
        inOutProperties{frameCrossing}=[inOutProperties{frameCrossing}; 0];
        association{frameCrossing}=[association{frameCrossing}; tIdx ];
        crossingEB3(frameCrossed).xCoord=[crossingEB3(frameCrossed).xCoord; [aTrack.x(crossIdx) 0.5]];
        crossingEB3(frameCrossed).yCoord=[crossingEB3(frameCrossed).yCoord; [aTrack.y(crossIdx) 0.5]];
        crossingEB3(frameCrossed).zCoord=[crossingEB3(frameCrossed).zCoord; [aTrack.z(crossIdx) 0.5]];
        inOutProperties{frameCrossed}=[inOutProperties{frameCrossed}; 1];
        association{frameCrossed}=[association{frameCrossed}; tIdx ];
    end
end
%%
cumulInterpCrossingEB3.xCoord=vertcat(cell2mat([arrayfun(@(x) x.xCoord,interpCrossingEB3,'unif',0)]));
cumulInterpCrossingEB3.yCoord=vertcat(cell2mat([arrayfun(@(x) x.yCoord,interpCrossingEB3,'unif',0)]));
cumulInterpCrossingEB3.zCoord=vertcat(cell2mat([arrayfun(@(x) x.zCoord,interpCrossingEB3,'unif',0)]));

% amiraWriteTracks('amiraTracks/tracks.am',tracks,'scales',[100,100,235],'edgeProp',{{'crossing',crossingTracks}});
% amiraWriteMovieInfo('amiraCrossingDetect/detectBeforeAfter.am',crossingEB3,'scales',[100,100,235],'prop',{{'inOut',inOutProperties},{'ID',association}});
% amiraWriteMovieInfo('amiraCrossCoordinate/crossing.am',interpCrossingEB3,'scales',[100,100,235]);
% amiraWriteMovieInfo('amiraCumulCrossCoordinate/cumulCrossing.am',cumulInterpCrossingEB3,'scales',[100,100,235]);

