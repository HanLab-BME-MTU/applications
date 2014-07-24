function results = statsuTrack(tracksFinal,varargin)
%STATSUTRACK performs analysis of STORM data using 'tracksFinal' of uTrack
%
% US, 01/10/2013
%

ip=inputParser;
ip.addRequired('tracksFinal',@isstruct);
ip.addOptional('rep',10,@isscalar);
ip.addOptional('tstep',100,@isscalar);
ip.addOptional('nFrames',3000,@isscalar);
ip.addOptional('tact',0.05,@isscalar);
ip.addOptional('texp',0.1,@isscalar);

ip.parse(tracksFinal,varargin{:});

tracksFinal=ip.Results.tracksFinal;
rep=ip.Results.rep;
tstep=ip.Results.tstep;
nFrames=ip.Results.nFrames;
tact=ip.Results.tact;
texp=ip.Results.texp;

nTracks=numel(tracksFinal);

tracks=repmat(...
    struct('trajectory',[],'driftCorrected',[],'centerOfMass',[],...
           'frameInfo',NaN(1,3),'fromUV',false,'gapClosed',false),nTracks,1);
       
trackLengths=zeros(nTracks,2);

for i=1:nTracks
       
    % determine start, end and length of trajectory
    seqOfEvents=tracksFinal(i).seqOfEvents;
    startFrame=min(seqOfEvents(:,1));
    endFrame=max(seqOfEvents(:,1));
    
    tracks(i).frameInfo(1)=startFrame;
    tracks(i).frameInfo(2)=endFrame;
    tracks(i).frameInfo(3)=endFrame-startFrame+1;
    
    trackLengths(i,:)=[endFrame-startFrame+1,startFrame];
    
    if( mod(startFrame,rep) == 1 )
        tracks(i).fromUV=true;
    end
    
    % positions, amplitudes, and frames of tracked particles
    p=tracksFinal(i).tracksCoordAmpCG;
    [nrows,ncols]=size(p);
    tmp=reshape(p,8,ncols/8)';
    tmp(:,end+1)=startFrame:endFrame;
    tracks(i).trajectory=tmp;
    
    idx=~isnan(tmp(:,1));
    tmp=tmp(idx,:);
    [nrows,ncols]=size(tmp);
    if( nrows > 1 )
        [wm,ws]=weightedStats(tmp(:,1:2),tmp(:,5:6));
        tracks(i).centerOfMass=[wm,ws,sqrt(sum(ws.^2))];
    else
        tracks(i).centerOfMass=[tmp(1:2), tmp(5:6), sqrt(sum(tmp(5:6).^2))];
    end
    
    % mark tracks with gap closing
    if( any(isnan(tmp(:))) )
        tracks(i).gapClosed=true;
    end
    
end

% analyse entire tracks
[histTotLen,xTotLen]=hist(trackLengths(:,1),1:max(trackLengths(:,1)));

% analyse tracks starting right after UV pulse
idx=vertcat(tracks.fromUV);
lengthUV=trackLengths(idx,:);
[histLenUV,xUV]=hist(lengthUV(:,1),1:max(lengthUV(:,1)));
% and those starting in the middle of an activation cycle
lengthNonUV=trackLengths(~idx,:);
[histLenNonUV,xNonUV]=hist(lengthNonUV(:,1),1:max(lengthNonUV(:,1)));

% save histograms as MATLAB figure

% calculate average track length dependent on starting frame
aveLengthUV=zeros(nFrames/tstep,3);
i=1;
k=1;
t=k*tstep;
tmp=[];
while( 1 )
    tmp=horzcat(tmp,lengthUV(i));
    
    if( lengthUV(i,2) > t )
        aveLengthUV(k,1)=t;
        aveLengthUV(k,2)=mean(tmp);
        aveLengthUV(k,3)=std(tmp);
        tmp=[];
        k=k+1;
        t=k*tstep;
    end
    
    i=i+1;
    if( i > size(lengthUV,1) )
        aveLengthUV(k,1)=t;
        aveLengthUV(k,2)=mean(tmp);
        aveLengthUV(k,3)=std(tmp);
        break;
    end
    
    if( t > nFrames )
        break;
    end
end

aveLengthNonUV=zeros(nFrames/tstep,3);
i=1;
k=1;
t=k*tstep;
tmp=[];
while( 1 )
    tmp=horzcat(tmp,lengthNonUV(i));
    
    if( lengthNonUV(i,2) > t )
        aveLengthNonUV(k,1)=t;
        aveLengthNonUV(k,2)=mean(tmp);
        aveLengthNonUV(k,3)=std(tmp);
        tmp=[];
        k=k+1;
        t=k*tstep;
    end
    
    i=i+1;
    if( i > size(lengthNonUV,1) )
        aveLengthNonUV(k,1)=t;
        aveLengthNonUV(k,2)=mean(tmp);
        aveLengthNonUV(k,3)=std(tmp);
        break;
    end
    
    if( t > nFrames )
        break;
    end
end

% estimated activation rates: UV-induced and spontaneous
konUV=size(lengthUV,1)/(tact*nFrames/rep);
konTH=size(lengthNonUV,1)/(texp*(nFrames-nFrames/rep));
% konTH=size(lengthNonUV,1)/(texp*nFrames);

% fit track length distributions to single exponential
%fitHistTot=fit(xTotLen(2:rep)',log(histTotLen(2:rep))','poly1');
%fitHistUV=fit(xUV(1:rep)',log(histLenUV(1:rep))','poly1');
%fitHistNonUV=fit(xNonUV(2:rep)',log(histLenNonUV(1:rep))','poly1');

% calculate frame-to-frame drift
shiftMatrix=NaN(nTracks,nFrames,2);
for k=1:nTracks
    shift=tracks(k).trajectory(:,1:2);
    startFrame=tracks(k).frameInfo(1);
    endFrame=tracks(k).frameInfo(2);
    tmp=[shift(1,1:2);shift(1:end-1,:)];
    shiftMatrix(k,startFrame:endFrame,:)=shift-tmp;
end

shift=squeeze(nanmean(shiftMatrix,1));
idx=isnan(shift(:,1));
shift(idx,:)=0.0;

% correct for drift
drift=cumsum(shift);
for k=1:nTracks
    startFrame=tracks(k).frameInfo(1);
    endFrame=tracks(k).frameInfo(2);
    t=tracks(k).trajectory;
    t(:,1:2)=t(:,1:2)-drift(startFrame:endFrame,:);
    tracks(k).driftCorrected=t;
end

% store results in struct 'results'
results.tracks=tracks;
results.shift=shift;
results.trackLengths=trackLengths;
results.histLenTot=[xTotLen',histTotLen'];
results.histLenUV=[xUV',histLenUV'];
results.histLenNonUV=[xNonUV',histLenNonUV'];
results.aveLengthUV=aveLengthUV;
results.aveLengthNonUV=aveLengthNonUV;
results.konUV=konUV;
results.konTH=konTH;
%results.fitHistTot=fitHistTot;
%results.fitHistUV=fitHistUV;
%results.fitHistNonUV=fitHistNonUV;