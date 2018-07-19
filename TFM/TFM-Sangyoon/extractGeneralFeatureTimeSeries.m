function meas=extractGeneralFeatureTimeSeries(curBcc)
% T=extractGeneralFeatureTimeSeries(curBcc) extracts.
% input:
%       curBcc: NxM matrix where N is the number of samples and M is the
%       maximum time limit
% output:
%       meas: features (NxL) matrix where L is the number of features
%       (currently L=8)
%       
% Sangyoon Han 2018.5.17

splineParam = 0.5;
%#1
maxIntensityNAs = nanmax(curBcc,[],2); %this should be high for group 2
endingSignal = zeros(size(maxIntensityNAs));
lifeTimeNAs = zeros(size(maxIntensityNAs));
timeToMaxInten = zeros(size(maxIntensityNAs));
ampSlopeNAs = zeros(size(maxIntensityNAs));
earlyAmpSlopeNAs = zeros(size(maxIntensityNAs));
lateAmpSlopeNAs = zeros(size(maxIntensityNAs));
%#5
meanIntensityNAs = nanmean(curBcc,2); %this should be high for group 2
for jj=1:size(curBcc,1)
    curData = curBcc(jj,:);
    curData(curData==0)=NaN;
    %2
    endingSignal(jj) = curData(find(~isnan(curData),1,'last')); %this should be high for group 2
    %#4
    lifeTimeNAs(jj) = sum(~isnan(curData)); %this should be low for group 6

    tRange = find(~isnan(curData));
    numNaNsInitial = find(~isnan(curData),1,'first')-1;
    d = curData(~isnan(curData));
    warning('off','SPLINES:CHCKXYWP:NaNs')
    sd_spline= csaps(tRange,d,splineParam);
    sd=ppval(sd_spline,tRange);
    % Find the maximum
    [~,curFrameMaxAmp]=nanmax(sd);
    timeToMaxInten(jj) = curFrameMaxAmp+numNaNsInitial;

    % 18-20. We have to use the amp slopes
    [~,ampSlopeNAs(jj)] = regression(tRange,d);
    earlyPeriod = 1:(min(10,length(d)));
    [~,earlyAmpSlopeNAs(jj)] = regression(tRange(earlyPeriod),d(earlyPeriod));
    latePeriod = (max(1,length(d)-9)):length(d);    
    [~,lateAmpSlopeNAs(jj)] = regression(tRange(latePeriod),d(latePeriod));
end

% meas=[maxIntensityNAs meanIntensityNAs endingSignal lifeTimeNAs ...
%     timeToMaxInten ampSlopeNAs earlyAmpSlopeNAs lateAmpSlopeNAs];
% meas=[timeToMaxInten lifeTimeNAs];% earlyAmpSlopeNAs lateAmpSlopeNAs];
meas=[timeToMaxInten];% lifeTimeNAs];% earlyAmpSlopeNAs lateAmpSlopeNAs];