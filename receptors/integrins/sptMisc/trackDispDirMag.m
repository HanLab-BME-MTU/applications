function trackDispDirMag(tracksFinal,diffAnalysisRes)

%get number of tracks
numTracks = length(tracksFinal);

%get number of segments making each track and the row of the first
%segment of each track if all track segments were put together
numSegPerTrack = zeros(numTracks,1);
for iTrack = 1 : numTracks
    numSegPerTrack(iTrack) = size(tracksFinal(iTrack).tracksCoordAmpCG,1);
end

%get total number of segments
numSegments = sum(numSegPerTrack);

%reserve memory for output parameters
[initCoord,finalCoord] = deal(NaN(numSegments,2));

%go over all compound tracks
for iTrack = 1 : numTracks
    
    %get number of segments in this compound track
    numSeg = numSegPerTrack(iTrack);
    
    %get current track's coordinates
    trackCoordCurrent = tracksFinal(iTrack).tracksCoordAmpCG;
    xCoord = trackCoordCurrent(:,1:8:end);
    yCoord = trackCoordCurrent(:,2:8:end);
    
    %get initial and final coordinates
    for jSeg = 1 : numSeg
        indxFirst = find(~isnan(xCoord(jSeg,:)),1,'first');
        initCoord(iTrack+jSeg-1,:) = [xCoord(jSeg,indxFirst) yCoord(jSeg,indxFirst)];
        indxLast = find(~isnan(xCoord(jSeg,:)),1,'last');
        finalCoord(iTrack+jSeg-1,:) = [xCoord(jSeg,indxLast) yCoord(jSeg,indxLast)];
    end
    
end

%calculate net displacement
dispNet = finalCoord - initCoord;

%get classifications from diffusion analysis results
trajClass = vertcat(diffAnalysisRes.classification);
trajClass = trajClass(:,2);

indxConf = find(trajClass==1);
indxBrown = find(trajClass==2);
indxDir = find(trajClass==3);
indxUndet = find(isnan(trajClass));

% imshow(ones(256))
hold on
plot(initCoord(indxConf,1),initCoord(indxConf,2),'c.')
quiver(initCoord(indxUndet,1),initCoord(indxUndet,2),dispNet(indxUndet,1),dispNet(indxUndet,2),'Color',[0.7 0.7 0.7])
quiver(initCoord(indxConf,1),initCoord(indxConf,2),dispNet(indxConf,1),dispNet(indxConf,2),'b')
quiver(initCoord(indxBrown,1),initCoord(indxBrown,2),dispNet(indxBrown,1),dispNet(indxBrown,2),'c')
quiver(initCoord(indxDir,1),initCoord(indxDir,2),dispNet(indxDir,1),dispNet(indxDir,2),'m')

%% ~~~ the end ~~~
