%SCRIPT testDiffStationarity generates many trajectories using mtGTPCapLDepK and tests for stationarity

%input to call mtGTPCapLDepK
initialState = struct('mtLength0',0.44,'capSize0',3,'unitConc',20);

modelParam = struct('minLength',0.2,'maxLength',1.5,'kTOnTFree',0.45,...
    'addAmpT',0.3,'addWidT',10,'addLenT',1.1,'kTOff',0.0,'kTOnD',0.2,...
    'kDOffFree',15,'addAmpD',10,'addWidD',10,'addLenD',0.6,'kHydrolysis',11);

maxNumSim = 1000;
totalTime = 149;
simTimeStep = 0.005;
timeEps = 0.5;
expTimeStep = 1;
aveInterval = 0.6;

saveTraj = struct('saveOrNot',0,'fileName',[]);

%get maxNumSim trajectories and average them to get specified sampling rate
mtLengthAve = zeros(floor((totalTime+expTimeStep)/expTimeStep),maxNumSim);
mtLengthSD = zeros(floor((totalTime+expTimeStep)/expTimeStep),maxNumSim);
for bigIter = 1:maxNumSim
    
    [mtLength,capSize,errFlag] = mtGTPCapLDepK(modelParam,initialState,...
        totalTime+expTimeStep,simTimeStep,timeEps,saveTraj);
    if errFlag
        return;
    end
    
    [mtLengthAve(:,bigIter),mtLengthSD(:,bigIter),errFlag] = averageMtTraj(...
        mtLength,simTimeStep,expTimeStep,aveInterval);
    if errFlag
        return;
    end
    
    mtLengthDiff(:,bigIter) = mtLengthAve(2:end,bigIter) - mtLengthAve(1:end-1,bigIter);
    
end

%save trajectories in file
save('trajectoriesEns','mtLengthDiff');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CHECK WHETHER PROCESS IS STATIONARY

%calculate the average difference and standard deviation at each time point 
mtDiffEnAve = mean(mtLengthDiff')'; 
mtDiffEnSD = std(mtLengthDiff')';

maxShift = 50;  %max tau in autocovarience function

%calculate the autocovarience function using ensemble averaging
enAveMat = repmat(mtDiffEnAve,1,maxNumSim);
autoCov = zeros(length(mtDiffEnAve)-2*maxShift,2*maxShift+1);
for shift = [-maxShift:1:maxShift]
    autoCov(:,maxShift+shift+1) = mean(...
        (mtLengthDiff(maxShift+1:end-maxShift,:)-enAveMat(maxShift+1:end-maxShift,:))'...
        .*...
        (mtLengthDiff(maxShift+shift+1:end-maxShift+shift,:)...
        -enAveMat(maxShift+shift+1:end-maxShift+shift,:))'...
        )';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CHECK WHETHER DISTRIBUTION IS GAUSSIAN

%get histogram approximating distribution of positions at each time step
nbins = 50;
[distribution,distVal] = hist(mtLengthDiff',nbins);
distribution = distribution*nbins/(distVal(end)+distVal(2)-2*distVal(1))/maxNumSim;

%get Gaussian distribution given the average position and its standard
%deviation at each time step
for i=1:length(mtDiffEnAve)
    gaussFunc(:,i) = exp(-(distVal-mtDiffEnAve(i)).^2/2/mtDiffEnSD(i)^2)...
        /mtDiffEnSD(i)/sqrt(2*pi);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%save analysis in file
save('dataEns','mtDiffEnAve','mtDiffEnSD','autoCov','distribution','distVal','gaussFunc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TEST FOR ERGODICITY

totalTime = maxNumSim-1;

%get trajectory and average it to get specified sampling rate
[mtLength,capSize,errFlag] = mtGTPCapLDepK(modelParam,initialState,...
    totalTime+expTimeStep,simTimeStep,timeEps,saveTraj);
if errFlag
    return;
end

[mtLengthAve,mtLengthSD,errFlag] = averageMtTraj(...
    mtLength,simTimeStep,expTimeStep,aveInterval);
if errFlag
    return;
end

mtLengthDiff = mtLengthAve(2:end) - mtLengthAve(1:end-1);

%save trajectory in file
save('trajectoriesT','mtLengthDiff');

%calculate the average length of MT and its standard deviation
meanL = mean(mtLengthDiff);
sdL = std(mtLengthDiff);

%calculate the autocovarience function using time averaging, to check
%whether system is ergodic
ergCov = xcov(mtLengthDiff,mtLengthDiff,maxShift,'biased'); 

%save analysis in file
save('dataT','meanL','sdL','ergCov');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
