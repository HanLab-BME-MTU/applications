%SCRIPT testModel tests the prediction ability of a model

totalTime = 1000;
simTimeStep = 0.005;
timeEps = 0.5;
expTimeStep = 1;
aveInterval = 0.6;

%input to call mtDynInstability
initialState = struct('mtLength0',0.44,'mtState0',1,'free',1,'unitConc',10);
modelParam = struct('minLength',-10000,'kOnElong',0.4,'kOffElong',5,...
    'kOnShrink',0.4,'kOffShrink',5);
%get trajectory using mtDynInstability
[mtLength,forcedRescue,errFlag] = mtDynInstability(modelParam,initialState,...
    totalTime+expTimeStep,simTimeStep,timeEps);
if errFlag
    return;
end

%average trajectory
[mtLengthAve,mtLengthSD,errFlag] = averageMtTraj(mtLength,simTimeStep,...
    expTimeStep,aveInterval);
if errFlag
    return;
end

lengthD1 = iddata(mtLengthAve,[],1);

[yh,fitV4(1)] = compare(lengthD1,model11,1);
[yh,fitV4(2)] = compare(lengthD1,model11,2);
[yh,fitV4(3)] = compare(lengthD1,model11,3);
[yh,fitV4(4)] = compare(lengthD1,model11,4);
[yh,fitV4(5)] = compare(lengthD1,model11,5);
fitV4 = fitV4';

