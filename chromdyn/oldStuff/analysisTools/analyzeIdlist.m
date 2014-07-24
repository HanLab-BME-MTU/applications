function anaDat=analyzeIdlist(idlist,dataProperties,options)
% ANALYZEIDLIST perform analysis of coordinates
% 
%
% SYNOPSIS anaDat=analyzeIdlist(idlist,dataProperties,options)
%
% INPUT idlist         : idlist or idlisttrack
%       dataProperties : structure
%       options        : select analysis (structure: options.process=[ 1 2 3 4 5]
%                                                     1 = tag separation, 2= relative translation,
%                                                     3 = relative difference, 4= fourier, 5= diffusion)
%                                          options.tags= [1 2]; select tags to be analyzed
%
% OUTPUT anaDat : analysisData structure

% c: 18/3/03	dT

%

if nargin<3 | isempty(options)
    options.process=[1 2 3 4 5];
    options.tags=[1 2];
end;
sTag=options.tags(1);
eTag=options.tags(2);
COLORLIST={'b','r-.','g:','k--','y'};
%init vars
gT=1;
sT=0;

acqTime=mean(dataProperties.frameTime,2);

%first nonempty entry
startT=1;
while(isempty(idlist(startT).linklist))
    startT=startT+1;
end;

numTags=size(idlist(startT).linklist,1);
anaDat(1).tags=options.tags;
anaDat(gT).coord=idlist(startT).linklist(:,9:11);
anaDat(gT).shift=[];
anaDat(gT).centroid=idlist(startT).centroid;
anaDat(gT).distMatrix=distMat(anaDat(gT).coord);
anaDat(gT).time=acqTime(startT);
anaDat(gT).sigSepChange=0;
% init temp data
[tempDat(1:length(idlist)).time]=deal(NaN);
[tempDat(1:length(idlist)).coord]=deal(zeros(numTags,3));

tempDat(startT).time=acqTime(startT);
tempDat(startT).coord=idlist(startT).linklist(:,9:11);

lt=startT;
lastSign=1;

%analyze it
for t=startT+1:length(idlist)
    if ~isempty(idlist(t).linklist)
        %temporary data for diffusion analysis
        tempDat(t).time=acqTime(t);
        tempDat(t).coord=idlist(t).linklist(:,9:11);
        
        anaDat(gT).shift=idlist(t).linklist(:,9:11)-anaDat(gT).coord;
        gT=gT+1;
        anaDat(gT).coord=idlist(t).linklist(:,9:11);
        anaDat(gT).centroid=idlist(t).centroid;
        anaDat(gT).distMatrix=distMat(anaDat(gT).coord);
        if squeeze(anaDat(gT-1).distMatrix(sTag,eTag))~=0 & squeeze(anaDat(gT).distMatrix(sTag,eTag))~=0 & isfield(idlist,'info')
            sig=testSigSeparationChange(idlist,dataProperties,[lt t],options.tags);
            anaDat(gT).sigSepChange=lastSign* sig;
            if sig~=0;
                lastSign=sig;
            end;
        else
             anaDat(gT).sigSepChange=0;
        end;
        anaDat(gT).time=acqTime(t);
        lt=t;
    end;
    
end;

%display results
coordShift=cat(3,anaDat.shift);
distMatrix=cat(3,anaDat.distMatrix);

if any(options.process==1)
    figure('Name','Separation & velocity')
    %plot distances
    separation=squeeze(distMatrix(sTag,eTag,:));
    sepNoFus=separation(find(separation));
    timeNoFus=[anaDat(find(separation)).time]';
    separation(find(separation==0))=nan;
    subplot(2,1,2);
    hold on;
    velocity=60*(sepNoFus(2:end)-sepNoFus(1:end-1))./(timeNoFus(2:end)-timeNoFus(1:end-1));
    plot((timeNoFus(2:end)+timeNoFus(1:end-1))/2,velocity,'b');
    plot((timeNoFus(2:end)+timeNoFus(1:end-1))/2,abs(velocity),'r');
    axis([0 max([anaDat.time]) min(velocity) max(abs(velocity))]);
    axes(gca);
    title(['spindle length velocity: ' idlist(1).stats.labelcolor{sTag} '-'  idlist(1).stats.labelcolor{eTag}]);
    xlabel('sec');
    ylabel('microns/minute');

    subplot(2,1,1);
    plot(timeNoFus,sepNoFus,'b:');
    hold on;
    plot([anaDat.time],separation,'r-');
    plot([anaDat(1).time anaDat(end).time],[mean(sepNoFus) mean(sepNoFus)],'k-');

    
    axis([0 max([anaDat.time]) 0 max(separation)]);
    axes(gca);
    title(['tag separation: ' idlist(1).stats.labelcolor{sTag} '-'  idlist(1).stats.labelcolor{eTag}]);
    xlabel('sec');
    ylabel('microns');
    
    t=[anaDat(find(isnan(separation)==0)).time]';
    X=[ones(size(t)) t];
    a=X\sepNoFus;
    plot(t,X*a);
    timeSig=[anaDat(find([anaDat.sigSepChange]==-1)-1).time];
    separation=squeeze(distMatrix(sTag,eTag,:));
    sepSig=separation(find([anaDat.sigSepChange]==-1)-1);
    plot(timeSig,sepSig,'ro');
end;

if any(options.process==2)
    figure('Name',['Relative Translation']);
    %number of coords
    for c=1:3
        %coords
        subplot(2,2,c);
        hold on;
        for i=2:numTags
            plot([anaDat(1:end-1).time],squeeze(coordShift(i,c,:)-coordShift(sTag,c,:)),COLORLIST{i-1});
        end;
        axes(gca);
        title(['x' num2str(c) '-Axis']);
        xlabel('sec');
        ylabel('microns');
    end
end;
%3D plot difference in translation
% 
% figure('Name',['Translation-diff']);
% trdiff=squeeze(coordShift(sTag,:,:)-coordShift(eTag,:,:));
% plot3(trdiff(1,:),trdiff(2,:),trdiff(3,:));

%plot difference in translation
if any(options.process==3)
    figure('Name',['Translation-diff']);
    trdiff=squeeze(coordShift(sTag,:,:)-coordShift(eTag,:,:))';
    subplot(2,1,1);
    plot([anaDat(1:end-1).time],sqrt(sum(trdiff.^2,2)));
    axes(gca);
    title(['length difference']);
    xlabel('sec');
    ylabel('microns');
    
    subplot(2,1,2);
    angle=squeeze(real(acos(dot(coordShift(sTag,:,:),coordShift(eTag,:,:),2)./sqrt(sum(coordShift(sTag,:,:).^2,2))./sqrt(sum(coordShift(eTag,:,:).^2,2)))));
    plot([anaDat(1:end-1).time],angle);
    
    axes(gca);
    title(['angle between translation vector']);
    xlabel('sec');
    ylabel('radian');
    
    %plot fft of some data
    
    meanacqTime=mean(diff(mean(dataProperties.frameTime,2)));
end;
if any(options.process==4)
    
    figure('Name',['FFT analysis']);
    coords=cat(3,anaDat.coord);
    coords1=squeeze(coords(sTag,:,:))';
    coords2=squeeze(coords(eTag,:,:))';
    fc1=fftn(coords1);nfc1=fc1.*conj(fc1);
    fc2=fftn(coords2);nfc2=fc2.*conj(fc2);
    separation=squeeze(distMatrix(1,2,:));
    fcs=fftn(separation);nfcs=fcs.*conj(fcs);
    
    freq=(2:length(nfc1)/2)/(length(nfc1)/2-1)*.5/meanacqTime;
    subplot(3,1,1);
    plot([freq' freq' freq'],nfc1(2:floor(end/2),:));
    axes(gca);
    title(['FFT of tag1 coords']);
    xlabel('frequency');
    
    freq=(2:length(nfc2)/2)/(length(nfc1)/2-1)*.5/meanacqTime;
    subplot(3,1,2);
    plot([freq' freq' freq'],nfc2(2:floor(end/2),:));
    axes(gca);
    title(['FFT of tag2 coords']);
    xlabel('frequency');
    
    freq=(2:length(nfcs)/2)/(length(nfc1)/2-1)*.5/meanacqTime;
    subplot(3,1,3);
    plot(freq,nfcs(2:floor(end/2)));
    axes(gca);
    title(['FFT of tag1-tag2 separation']);
    xlabel('frequency');
end;
%diffusion analysis
for s=1:floor(length(anaDat)/2)
    deltaPos=cat(3,tempDat(1+s:end).coord)-cat(3,tempDat(1:end-s).coord);
    deltaPosTags=deltaPos(eTag,:,:)-deltaPos(sTag,:,:);
    deltaT=[tempDat(1+s:end).time] -[tempDat(1:end-s).time];
    validIdx=find(isnan(deltaT)==0);
    dL(:,s)=mean((sum(deltaPos(:,:,validIdx).^2,2)),3);
    dLT(s)=mean((sum(deltaPosTags(:,:,validIdx).^2,2)),3);
end
%plot diffusion
dT=dataProperties.timeLapse:dataProperties.timeLapse:dataProperties.timeLapse*(size(dL,2));
if any(options.process==5)
    figure('Name',['Diffusion']);
    hold on;
    for tg=1:numTags
        h(tg)=plot(dT,dL(tg,:),COLORLIST{tg});
    end;
    % add relative pos
    h(tg+1)=plot(dT,dLT,'k-');
    legend(h,{idlist(1).stats.labelcolor{1:numTags} [idlist(1).stats.labelcolor{sTag} '-'  idlist(1).stats.labelcolor{eTag}]});
    %axes(gca);
    title(['Diffusion']);
    xlabel('{\DeltaT [sec]}');
    ylabel('r^2 [micron]');
end;