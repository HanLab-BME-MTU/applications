%
% an approach to count single molecules
% 
% supplied input variables on workspace
%   -> tracksFinal as output from uTrack package
%   
% US, 2013/01/30
%

% some basic definitions that are likely to be changed
%   tExp  ->  exposure time in seconds
% imSize  ->  image size in pixel
%    rep  ->  number of frames in activation cycle
tExp=0.05;
imSize=[256,256];
rep=10;

nTracks=numel(tracksFinal);

%
% new struct 'track' with the following members:
%
%      coord  ->  coordinates and corresponding errors
%        com  ->  center of mass
%   timeInfo  ->  start, end, duration in frames
%  gapClosed  ->  gap closing used?
%     fromUV  ->  induced by UV activation?
%
track=repmat(...
    struct('coord',[],'com',[],'amp',[],'timeInfo',[],'gapClosed',false,'fromUV',false),nTracks,1);

% index of UV induced tracks
idxUV=false(nTracks,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                             
%%%  bring tracksFinal in human-readable form
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:nTracks
    
    startF=min( tracksFinal(i).seqOfEvents(:,1) );
    stopF=max( tracksFinal(i).seqOfEvents(:,1) );
    
    track(i).timeInfo(1)=startF;
    track(i).timeInfo(2)=stopF;
    track(i).timeInfo(3)=stopF-startF+1;
    
    % mark tracks right after UV activation
    if mod(startF,rep) == 1
        track(i).fromUV=true;
        idxUV(i)=true;
    end
    
    % coordinates, amplitudes and corresponding errors
    p=tracksFinal(i).tracksCoordAmpCG;
    [nrows,ncols]=size(p);
    tmp=reshape(p,8,ncols/8)';
    tmp(:,end+1)=[startF:stopF];
    track(i).coord=tmp(:,[1,2,5,6]);
    track(i).amp=tmp(:,[4,8,end]);
    % center of mass calculated as a weigthed mean
    idx=~isnan(tmp(:,1));
    tmp=tmp(idx,:);
    [nrows,ncols]=size(tmp);
    if( nrows > 1 )
        [wm,ws]=weightedStats(tmp(:,1:2),tmp(:,5:6),'s');
        track(i).com=[wm,ws,sqrt(sum(ws.^2))];
    else
        track(i).com=[tmp(1:2), tmp(5:6), sqrt(sum(tmp(5:6).^2))];
    end
    % mark klusters with gap closing
    if any( isnan(tmp(:)) )
        track(i).gapClosed=true;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  determine rates of the underlying photo-kinetics
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% start with distribution of ON times
timeON=vertcat(track.timeInfo);
timeON=timeON(:,3)*tExp;
bins=0:tExp:max(timeON);
[yON,xON]=hist(timeON,bins);

% create and fit CDF of ON times: 1-exp(-kc*x)
% kc is the sum of the OFF and bleaching rate
area=tExp*sum(yON);
yONcdf=cumsum(yON)*tExp/area;

ft=fittype(@(kc,x) 1.0-exp(-kc*x));
pini=( yONcdf(2)-yONcdf(1) )/( xON(2)-xON(1) );
tonFit=fit(xON(3:end)',yONcdf(3:end)',ft,'StartPoint',pini);

%% continue with distribution of blinks

% group center of mass of tracks via mean-shift clustering
allCom=vertcat(track.com);
[clusterInfo,clusterMap]=MeanShiftClustering(allCom(:,1:2),0.25,'kernel','flat');

for i=1:numel(clusterInfo)
    numPts=clusterInfo(i).numPoints;
    tmp=[];
    tmp2=[];
    for k=1:numPts
        ptID=clusterInfo(i).ptIdData(k);
        tmp=vertcat(tmp,[track(ptID).coord track(ptID).amp]);
        tmp2=vertcat(tmp2,track(ptID).timeInfo);
    end
    clusterInfo(i).data=tmp;
    clusterInfo(i).time=tmp2;
    
    if size(tmp,1) > 1
        [wm,ws]=weightedStats(tmp(:,1:2),tmp(:,3:4),'s');
        clusterInfo(i).ptClusterCenter=[wm,ws,sqrt(sum(ws.^2))];
    else
        clusterInfo(i).ptClusterCenter=[tmp(1:4),sqrt(sum(tmp(3:4).^2))];
    end
end

%% number of blinks, fit to CDF of geometric distribution
nBlinks=vertcat(clusterInfo.numPoints);
nBlinks=nBlinks-1;
[yBl,xBl]=hist(nBlinks,0:max(nBlinks));

area=sum(yBl);
yBlcdf=cumsum(yBl)/area;

% parameter of geometric distribution q = koff/(koff+kbleach)
ft=fittype('1.0-q^(x+1)');
pini=0.5;
blinkFit=fit(xBl',yBlcdf',ft,'StartPoint',pini);

%% calculate koff and kbleach from q and kc
q=blinkFit.q;
kc=tonFit.kc;

koff=q*kc;
kbleach=kc-koff;

% error propagation
tmp1=confint(blinkFit);
q_err=q-tmp1(1);
tmp2=confint(tonFit);
kc_err=kc-tmp2(1);
koff_err=sqrt( (kc*q_err)^2 + (q*kc_err)^2 );
kbleach_err=sqrt( koff_err^2 + kc_err^2 );


%% continue with distribution of OFF times
toff=[];
for k=1:numel(clusterInfo)
    timeInfo=clusterInfo(k).time;
    toff=vertcat(toff,timeInfo(2:end,1)-timeInfo(1:end-1,2));
    toff=vertcat(toff,timeInfo(1,1));
end

toff=toff*tExp;
[yOFF,xOFF]=hist(toff,0:tExp:max(toff));
yOFFcdf=cumsum(yOFF)/sum(yOFF);

ft=fittype(@(k1,x) 1.0-exp(-k1*x));
pini=( yOFFcdf(2)-yOFFcdf(1) )/( xOFF(2)-xOFF(1) );
toffFit=fit(xOFF',yOFFcdf',ft,'StartPoint',pini);

%% temporal dissection of clusters found via mean-shift


%% clear unused variables
clear startF stopF idx tmp i;
clear p nrows ncols nTracks;
clear tmp1 tmp2;