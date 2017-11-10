function [ track ] = reformTracksFinal(tracksFinal, driftLen)
%reformTracksFinal takes the output of Khulouds tracker and reformats it
%into a managable format.
%
% Inputs : 
%          tracksFinal, the output of the tracking program
%          driftLen, the lenght of a track inorder to considered a drift
%                     marker
%
%
% Outputs :
%           tracks, the refromated final product
%
% new struct 'track' with the following members:
%
%      coord  ->  coordinates and corresponding errors [x,y,dx,dy]
%        com  ->  center of mass [amp,d_amp, ?]
%   timeInfo  ->  [start, end, duration in frames, number of UvActivations]
%  gapClosed  ->  gap closing used?
%     fromUV  ->  induced by UV activation?
%    isDrift  ->  is track a drift marker (lenght > driftLen)
%
%
%
% written by Jeffrey Werbin adapted from CountVonCount 2013/07/15

nTracks=numel(tracksFinal);

maxGap = 10;
rep = 10;

% new struct 'track' with the following members:
%
%      coord  ->  coordinates and corresponding errors
%        com  ->  center of mass 
%   timeInfo  ->  start, end, duration in frames
%  gapClosed  ->  gap closing used?
%     fromUV  ->  induced by UV activation?
%     isDrift ->  is this a drift marker
%
track=repmat(...
    struct('coord',[],'com',[],'amp',[],'timeInfo',[],'gapClosed',false,'fromUV',false,'isDrift',false,'num',[]),nTracks,1);

% index of UV induced tracks
idxUV=false(nTracks,1);

% an index of which tracks span the whole movie
isDrift = false(nTracks,1);

% an index of which tracks have more than one point
isMany = false(nTracks,1);

%
isearly = false(nTracks,1);

%catches crazy bugs
isweird = false(nTracks,1);



for i=1:nTracks
    
    startF=min( tracksFinal(i).seqOfEvents(:,1) );
    stopF=max( tracksFinal(i).seqOfEvents(:,1) );
    
%     if (startF > MD.nFrames_) | (stopF > MD.nFrames_)
%         iswierd(i) = true;
%         continue;
%     end
    
    track(i).timeInfo(1)=startF;
    track(i).timeInfo(2)=stopF;
    track(i).timeInfo(3)=stopF-startF+1;
    track(i).timeInfo(4)=floor(startF/10)+1; %number of uv activations it has been exposed to
        
    if startF < 11
        isearly(i)=true;
    end
    
    % mark tracks right after UV activation
    if mod(startF,rep) == 1
        track(i).fromUV=true;
        idxUV(i)=true;
    end
    
    %mark tracks if they span the whole movie
    if track(i).timeInfo(3) >= driftLen
        track(i).isDrift=true;
    end
    
    % coordinates, amplitudes and corresponding errors
    p=tracksFinal(i).tracksCoordAmpCG;
    [nrows,ncols]=size(p);
    tmp=reshape(p,8,ncols/8)';
    tmp(:,end+1)=[startF:stopF];
    track(i).coord=tmp(:,[1,2,5,6]);
    track(i).amp=tmp(:,[4,8,end]);
    track(i).num = numel(tmp(:,1));
    % center of mass calculated as a weigthed mean
    idx=~isnan(tmp(:,1));
    tmp=tmp(idx,:);
    [nrows,ncols]=size(tmp);
    if( nrows > 1 )
        [wm,ws]=weightedStats(tmp(:,1:2),tmp(:,5:6),'s');
        track(i).com=[wm,ws,sqrt(sum(ws.^2))];
        isMany(i) = true;
        
        %approximates refited amplitudes after merging all time points and
        %refitting
        track(i).TotAmp = sum(track(i).amp(~isnan(track(i).amp(:,1)),1))/sqrt(sum(~isnan(track(i).amp(:,1))));
        
    else
        track(i).com=[tmp(1:2), tmp(5:6), sqrt(sum(tmp(5:6).^2))];
        track(i).TotAmp = track(i).amp(1,1);
    end
    
    % mark klusters with gap closing
    if any( isnan(tmp(:)) )
        track(i).gapClosed=true;
    end
end


end

