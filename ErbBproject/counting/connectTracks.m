function tracksCombined=connectTracks(tracks,varargin)
%CONNECTTRACKS links tracks from single molecules in space and time
%   
%
% US, 2013/21/01
%

% check input
ip=inputParser;

ip.addRequired('tracks',@isstruct);
ip.addOptional('ton',2.5,@isnumeric);
ip.addOptional('toff',10.0,@isnumeric);

ip.parse(tracks,varargin{:});

tracks=ip.Results.tracks;
ton=ip.Results.ton;
toff=ip.Results.toff;

nTracks=numel(tracks);

% output struct array
tracksCombined=struct();

% build KD tree with center of mass of each track
%ptlist=vertcat(tracks.centerOfMass);
%kd=KDTree(ptlist(:,1:2));

% all tracks starting inf frame 1
startEnd=vertcat(tracks.frameInfo);
idxFirst=startEnd(:,1) == 1;

% init first entries of tracksCombined with tracks starting in frame 1
k=1;
nCombTracks=0;
while( idxFirst(k) )
    path=tracks(k).trajectory;
    tracksCombined(k).path=path;
    tracksCombined(k).info=tracks(k).frameInfo;
    
    tracksCombined(k).centerOfMass=tracks(k).centerOfMass;
    nCombTracks=nCombTracks+1;
    tracksCombined(k).id=nCombTracks;
    tracksCombined(k).howMany=1;
    
    k=k+1;
end

% distribute remaining tracks
while( k <= nTracks )
    
    % information of trajectory to be distributed
    currentPath=tracks(k).trajectory;
    %idx=~isnan(currentPath(:,1));
    currentInterval=currentPath(:,end);
    
    currentMean=tracks(k).centerOfMass(1:2);
    %currentStd=tracks(k).centerOfMass(3:4);
    currentDelta=tracks(k).centerOfMass(5);
    currentDeltaSq=currentDelta*currentDelta;
    
    % check spatial proximity
    allMeans=vertcat(tracksCombined.centerOfMass);
    delta=allMeans(:,end);
    deltaSq=delta.*delta;
    allMeans=allMeans(:,1:2);
    squaredDistances=distance(allMeans',currentMean');
    idxDist=squaredDistances < deltaSq+currentDeltaSq;
    closeTrackID=vertcat(tracksCombined(idxDist).id);
    
    % current track is remote from all present tracks -> create new
    if( sum(idxDist) == 0 )
        nCombTracks=nCombTracks+1;
        tracksCombined(nCombTracks).path=currentPath;
        tracksCombined(nCombTracks).centerOfMass=tracks(k).centerOfMass;
        tracksCombined(nCombTracks).info=tracks(k).frameInfo;
        tracksCombined(nCombTracks).id=nCombTracks;
        tracksCombined(nCombTracks).howMany=1;
    else
        % check temporal overlap of current and present tracks
        allIntervals=repmat(struct('frames',[]),nCombTracks,1);
        for l=1:nCombTracks
            tmp=tracksCombined(l).info;
            for m=1:size(tmp,1)
                allIntervals(l).frames=...
                    horzcat(allIntervals(l).frames,[tmp(m,1):tmp(m,2)]);
            end
        end
        idxFrame=tempOverlap(currentInterval,allIntervals);
        clear allIntervals;
        
        % current track overlaps with each present track -> create new
        if( sum(idxFrame) == 0 )
            nCombTracks=nCombTracks+1;
            tracksCombined(nCombTracks).path=currentPath;
            tracksCombined(nCombTracks).centerOfMass=tracks(k).centerOfMass;
            tracksCombined(nCombTracks).info=tracks(k).frameInfo;
            tracksCombined(nCombTracks).id=nCombTracks;
            tracksCombined(nCombTracks).howMany=1;
       
        % current track overlpas with all but one track -> merge these
        elseif( sum(idxFrame) == 1 )
            % get track ID
            trackID=closeTrackID(idxFrame);
            
            tmp=tracksCombined(trackID).path;
            tmp=[tmp;currentPath];
            tracksCombined(trackID).path=tmp;
            
            idx=~isnan(tmp(:,1));
            tmp=tmp(idx,:);
            [wm,ws]=weightedStats(tmp(:,1:2),tmp(:,5:6));
            tracksCombined(trackID).centerOfMass=[wm,ws,sqrt(sum(ws.^2))];
            
            tmp=tracksCombined(trackID).info;
            tmp=[tmp; tracks(k).frameInfo];
            tracksCombined(trackID).info=tmp;     
            
            tracksCombined(trackID).howMany=...
                tracksCombined(trackID).howMany+1;
            
        % current track overlaps with more than one present track -> merge with closest    
        elseif( sum(idxFrame) > 1 )
            tmp=squaredDistances;
            tmp(~idxFrame)=NaN;
            [~,idx]=min(tmp);
            
            trackID=idx;%closeTrackID(idx);
            
            tmp=tracksCombined(trackID).path;
            tmp=[tmp;currentPath];
            tracksCombined(trackID).path=tmp;
            
            idx=~isnan(tmp(:,1));
            tmp=tmp(idx,:);
            [wm,ws]=weightedStats(tmp(:,1:2),tmp(:,5:6));
            tracksCombined(trackID).centerOfMass=[wm,ws,sqrt(sum(ws.^2))];
            
            tmp=tracksCombined(trackID).info;
            tmp=[tmp; tracks(k).frameInfo];
            tracksCombined(trackID).info=tmp;
            tracksCombined(trackID).howMany=...
                tracksCombined(trackID).howMany+1;
        end
    end
    
    k=k+1;
end


end

function idxFrame=tempOverlap(fseq,allSeq)

idxFrame=false(size(allSeq,1),1);

for n=1:size(allSeq,1)
    idxFrame(n)=isempty(intersect(fseq,allSeq(n).frames));
end

end