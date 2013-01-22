function tracksCombined=connectTracks(tracks,varargin)
%CONNECTTRACKS links single tracks from single molecules in space and time
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

% output struct array
tracksCombined=struct();

startEnd=vertcat(tracks.info);
idxFirst=startEnd(:,1) == 1;

% init first entries of tracksCombined with tracks starting in frame 1
k=1;
nCombTracks=0;
while( idxFirst(k) )
    path=tracks(k).trajectory;
    tracksCombined(k).path=path;
    tracksCombined(k).info=tracks(k).info;
    
    % calculate weigthed mean of trajectory
    idx=~isnan(path(:,1));
    [wm,ws]=weigthedStats(path(idx,1:2),path(idx,5:6));
    
    tracksCombined(k).center=[wm,ws,sqrt(sum(ws.^2))];
    nCombTracks=nCombTracks+1;
    
    k=k+1;
end

% distribute remaining tracks
while( true )
    
    % information of trajectory to be distributed
    currentPath=tracks(k).trajectory;
    idx=~isnan(currentPath(:,1));
    currentInterval=currentPath(:,end);
    
    [currentMean,currentStd]=weigthedStats(currentPath(idx,1:2),currentPath(idx,5:6));
    currentDelta=sum(currentStd.^2);
    
    % check temporal overlap
    allIntervals=vertcat(tracksCombined.info);
    idxInterval=tempOverlap(currentInterval,allIntervals);
    
    % check spatial proximity
    allMeans=vertcat(tracksCombined.center);
    delta=allMeans(:,end);
    delta=delta*delta;
    allMeans=allMeans(:,1:2);
    squaredDistances=distance(allMeans',currentMean');
    idxDist=squaredDistances < delta+currentDelta;
    
    % path holds current trajectory
    %path=tracks(k).trajectory;
    %idx=~isnan(path(:,1));
    %[wm,ws]=weigthedStats(path(idx,1:2),path(idx,5:6));
    %delta=sqrt(sum(ws.^2));
    % centers obtained so far
    %centers=vertcat(tracksCombined.center);
    %delta=delta+centers(:,end);
    %centers=centers(:,1:2);
    % squared distance between current centers and present mean
    %distM=distance(centers',wm');
    %idxDist=distM < delta*delta;
    
    
    % temporal overlap of trajectories
    %startEnd=path(:,end);
    %allStartEnd=vertcat(tracksCombined.info);
    %idxFrame=tempOverlap(startEnd,allStartEnd);
    
    %idx=idxDist & idxFrame;
    
%     if( sum(idx) == 0 )
%         tracksCombined(k).path=path;
%         tracksCombined(k).center=[wm,ws,sqrt(sum(ws.^2))];
%         nCombTracks=nCombTracks+1;
%         k=k+1;
%         continue;
%     elseif( sum(idx) == 1 )
%         l=find(idx == 1);
%         tmp=tracksCombined(l);
%         tmp=vertcat(tmp,path);
%         tracksCombined(l).path=tmp;
%         [wm,ws]=weightedStats(tmp(:,1:2),tmp(:,5:6));
%         tracksCombined(l).center=[wm,ws,sqrt(ws.^2)];
%     end
    
    
end


end

function idxFrame=tempOverlap(fseq,allSeq)

idxFrame=false(size(allSeq,1),1);

for n=1:size(allSeq,1)
    first=allSeq(n,1);
    last=allSeq(n,2);
    idxFrame(n)=isempty(intersect(fseq,first:last));
end

end