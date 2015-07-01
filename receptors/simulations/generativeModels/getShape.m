function [simDuration, corDuration, actualDuration] = getShape(winSize,lagT,locError)

kOff = 1;
N = 100;
dT = .01;

%get the simulated correlations of N sets of trajectors;
[simCorr,corr,actualDuration] = troubleshootSimulateProb(kOff,N,winSize,lagT,locError);

[simDuration,~] = fitCorr(simCorr,winSize,lagT,N);
[corDuration,~] = fitCorr(corr,winSize,lagT,N);
[~,numFrames] = size(corr);

kOffSim = getKOff(simDuration(simDuration>0),numFrames*.9)/dT;
kOffCor = getKOff(corDuration(corDuration>0),numFrames*.9)/dT;
kOffActual = getKOff(actualDuration,max(actualDuration)+2)/dT;
figure
maxVal = max(max(simDuration),max(actualDuration));
plot(actualDuration,simDuration,'r.',[0 maxVal],[0 maxVal],'k-')
xlabel('Actual bound time')
ylabel('Measured bound time')
title(sprintf('Simulation: kOff=%g, kOff-Actual=%g, Error=%gnm, winSize=%d, lag=%d',kOffSim,kOffActual,locError*1000,winSize,lagT))

figure
maxVal = max(max(corDuration),max(actualDuration));
plot(actualDuration,corDuration,'r.',[0 maxVal],[0 maxVal],'k-')
xlabel('Actual bound time')
ylabel('Measured bound time')
title(sprintf('Corr Coeff: kOff=%g, kOff-Actual=%g, Error=%gnm, winSize=%d, lag=%d',kOffCor,kOffActual,locError*1000,winSize,lagT))



function [duration,startT] = fitCorr(corr,winSize,lagT,N)
[~,numFrames] = size(corr);

winLag = winSize+lagT;

startT = zeros(N,1);
duration = zeros(N,1);
side = winLag - 2;

for n = 1:N
    fprintf('\n%d',n);
    minError = realmax('double');
    numGoodPoints = sum(corr(n,:)>0.2);
    approxH = sum(corr(n,:).*(corr(n,:)>0.2))/numGoodPoints;
    for h = (approxH*.9):0.01:(approxH*1.1)
        slope = h/side;
        for width = 1:min(round(numGoodPoints*1.1),numFrames)
            for offset = -1:numFrames
                trapezoid = zeros(1,numFrames);
                for j = max(offset,2):min((offset+2*side+width-1),numFrames)
                    if j <= (offset+side-1)
                        trapezoid(j) = trapezoid(j-1) + slope;
                    elseif j <= (offset+side+width-1)
                        trapezoid(j) = h;
                    else
                        trapezoid(j) = trapezoid(j-1) - slope;
                    end
                end
                error = sum((trapezoid-squeeze(corr(n,:))).^2);
                if error < minError
                    minError = error;
                    startT(n) = offset+2;
                    if (width+winLag<numFrames)
                        duration(n) = width + winLag;
                    else
                        duration(n) = numFrames+1;
                    end
                end
            end
        end
    end
end


function bestK = getKOff(data,maxVal)

data = sort(data);
cdf = zeros(sum(data<maxVal),2);

if length(cdf) >= length(data)/2
    cdf(:,1) = data(1:length(cdf));
    cdf(:,2) = (1:length(cdf))/length(data);
    cdf = [[0 0];cdf];
    minError = realmax('double');
    bestK = 0;
    for k = .01/max(data):.01/max(data):1
        fitCDF = 1 - exp(-k*(cdf(:,1)));
        error = sum((fitCDF-cdf(:,2)).^2);
        if error<minError
            bestK = k;
            minError = error;
        end
    end
else
    bestK = log(2)/median(data);
end