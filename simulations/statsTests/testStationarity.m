%SCRIPT testStationarity generates many trajectories using mtGTPCapLDepK and tests for stationarity

%input to call mtGTPCapLDepK
initialState = struct('mtLength0',0.44,'capSize0',3,'unitConc',20);

modelParam = struct('minLength',0.2,'maxLength',1.5,'kTOnTFree',0.45,...
    'addAmpT',0.3,'addWidT',10,'addLenT',1.1,'kTOff',0.0,'kTOnD',0.2,...
    'kDOffFree',15,'addAmpD',10,'addWidD',10,'addLenD',0.6,'kHydrolysis',11);

maxNumSim = 10000;
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
    
end

%save trajectories in file
%save('trajectoriesEns','mtLengthAve','mtLengthSD');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CHECK WHETHER PROCESS IS STATIONARY

%calculate the average length and standard deviation at each time point 
mtLengthEnAve = mean(mtLengthAve')'; 
mtLengthEnSD = std(mtLengthAve')';

maxShift = 50;  %max tau in autocovarience function

%calculate the autocovarience function using ensemble averaging
enAveMat = repmat(mtLengthEnAve,1,maxNumSim);
autoCov = zeros(length(mtLengthEnAve)-2*maxShift,2*maxShift+1);
for shift = [-maxShift:1:maxShift]
    autoCov(:,maxShift+shift+1) = mean(...
        (mtLengthAve(maxShift+1:end-maxShift,:)-enAveMat(maxShift+1:end-maxShift,:))'...
        .*...
        (mtLengthAve(maxShift+shift+1:end-maxShift+shift,:)...
        -enAveMat(maxShift+shift+1:end-maxShift+shift,:))'...
        )';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CHECK WHETHER DISTRIBUTION IS GAUSSIAN

%get histogram approximating distribution of positions at each time step
nbins = 50;
[distribution,distVal] = hist(mtLengthAve',nbins);
distribution = distribution*nbins/(distVal(end)+distVal(2)-2*distVal(1))/maxNumSim;

%get Gaussian distribution given the average position and its standard
%deviation at each time step
for i=1:length(mtLengthEnAve)
    gaussFunc(:,i) = exp(-(distVal-mtLengthEnAve(i)).^2/2/mtLengthEnSD(i)^2)...
        /mtLengthEnSD(i)/sqrt(2*pi);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%save analysis in file
%save('dataEns','mtLengthEnAve','mtLengthEnSD','autoCov','distribution','distVal','gaussFunc');

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

%save trajectory in file
%save('trajectoriesT','mtLengthAve','mtLengthSD');

%calculate the average length of MT and its standard deviation
meanL = mean(mtLengthAve);
sdL = std(mtLengthAve);

%calculate the autocovarience function using time averaging, to check
%whether system is ergodic
ergCov = xcov(mtLengthAve,mtLengthAve,maxShift,'biased'); 

%save analysis in file
%save('dataT','meanL','sdL','ergCov');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
