function [trajAve,trajSD,autoCov,ergCov,distribution,distVal,gaussFunc,H]...
    = testSNE(enSize,trajEn,trajT)
%TESTSNE tests a trajectory for weak stationarity, normality and ergodicity
%
%SYNOPSIS [trajAve,trajSD,autoCov,ergCov,distribution,distVal,gaussFunc,H]...
%    = testSNE(enSize,trajEn,trajT)
%
%INPUT  enSize  : Number of trajectories in ensemble, which equals length
%                 of trajectory "trajT".
%       trajEn  : Ensemble of trajectories used for ensemble averaging:
%                 # of columns = trajectory length, # of rows = enSize.
%       trajT   : Column vector (of size enSize) of 1 long trajectory, 
%                 used for time averaging.
%
%OUTPUT trajAve      :
%       trajSD       :
%       autoCov      :
%       ergCov       :
%       distribution :
%       distVal      :
%       gaussFunc    :
%       H            :
%
%Khuloud Jaqaman, February 2004

%calculate ensemble average length and standard deviation at each time point 
trajEnAve = mean(trajEn')'; 
trajEnSD = std(trajEn')';

%calculate the time average and standard deviation of trajectory
meanL = mean(trajT);
sdL = std(trajT);

maxShift = 50;  %max tau in autocovarience function

%calculate the autocovarience function using ensemble averaging
enAveMat = repmat(trajEnAve,1,enSize);
autoCov = zeros(length(trajEnAve)-2*maxShift,2*maxShift+1);
for shift = [-maxShift:1:maxShift]
    autoCov(:,maxShift+shift+1) = mean(...
        (trajEn(maxShift+1:end-maxShift,:)-enAveMat(maxShift+1:end-maxShift,:))'...
        .*(trajEn(maxShift+shift+1:end-maxShift+shift,:)...
        -enAveMat(maxShift+shift+1:end-maxShift+shift,:))')';
end

%calculate the autocovarience function using time averaging
ergCov = xcov(trajT,trajT,maxShift,'biased'); 

%get histogram approximating distribution of positions at each time step
nbins = 50;
[distribution,distVal] = hist(trajEn',nbins);
distribution = distribution*nbins/(distVal(end)+distVal(2)-2*distVal(1))/enSize;

%get Gaussian distribution given the average position and its standard
%deviation at each time step
for i=1:length(trajEnAve)
    gaussFunc(:,i) = exp(-(distVal-trajEnAve(i)).^2/2/trajEnSD(i)^2)...
        /trajEnSD(i)/sqrt(2*pi);
end

%test for normality
for i=1:length(trajEnAve)
    H(i) = lillietest(trajEn(i,:));
end

%put together ensemble average and time average
trajAve = [trajEnAve; meanL];
trajSD = [trajEnSD; sdL];
