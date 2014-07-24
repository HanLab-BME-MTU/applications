function [inputDataEntry,errFlag] = getAnaDat(mtLength,mtLengthErr,...
    simNum,timeStep,aveInterval)
%GETANADAT writes MT trajectory data in format required for calcMTDynamics. 

%SYNOPSIS [inputDataEntry,errFlag] = getAnaDat(mtLength,mtLengthErr,...
%    simNum,timeStep,aveInterval)
%
%INPUT mtLength    : MT length at every time point. 
%      mtLengthErr : Error in MT length at every time point.
%      simNum      : Simulation number.
%      timeStep    : Step between consecutive time points. 
%      aveInterval : Interval about time point represented by data.
%
%OUPUT inputDataEntry : Structure providing data about trajectory and simulation,
%                       in the format required by calcMTDynamics.
%      errFlag        : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman 10/03

errFlag = 0;

vecLength = length(mtLength);

simData(1:vecLength)=struct('info',[],'coord',[],'centroid',[],...
    'distanceMatrix',[],'distanceVectorMatrixN',[],'displacementVectorN',[],...
    'displacement',[],'time',[],'timepoint',[],'sigmaTime',[],'stats',[]);

%general information
simData(1).info.name = sprintf('simMtDynamicsRun%i',simNum);
simData(1).info.idlisttype = 'synthetic';
simData(1).info.nTags = 2;
simData(1).info.labelColor = [{'spb'};{'mtTip'}];
simData(1).info.created = nowString;
simData(1).info.dataProperties = [];

warning off MATLAB:divideByZero;

%data
for i=1:vecLength-1
    
    lengthI = mtLength(i); %length of MT at i (micrometers)
    simData(i).coord = [0 0 0; lengthI 0 0]; %position of SPB and MT-tip at i 
    %(micrometers)
    simData(i).centroid = [lengthI/2 0 0]; %position of center of geometry (micrometers)
    [simData(i).distanceMatrix, distanceVectorMatrix] = ... %matrix of distances between  
        distMat(simData(i).coord); %SPB&SPB, SPB&MT-tip, MT-tip&SPB, MT-tip&MT-tip, and
    %matrices of vectors pointing from one to the other
    simData(i).distanceVectorMatrixN = ... %normalized "distanceVectorMatrix"
        distanceVectorMatrix./repmat(simData(i).distanceMatrix,[1 1 3]);
    
    simData(i).displacementVectorN = [0 0 0; 1 0 0];  %of displacement
    simData(i).displacement = [0;mtLength(i+1)-lengthI]; %magnitude of displacement
    %between i and i+1 (micrometers)
    
    simData(i).time = (i-1)*timeStep; %time at iteration i (seconds)
    simData(i).timepoint = i; %iteration number, i.e.
    simData(i).sigmaTime = aveInterval/2; %interval of averaging (seconds) over 2
    
    simData(i).stats.qMatrix = diag([0.000064 0 0 mtLengthErr(i)^2 0 0 ...
            (0.000064+mtLengthErr(i)^2)/8 0 0]); %variance in position of SPB, MT-tip 
    %and centroid (micrometers^2)
    simData(i).stats.noise = [1 1];
    
end

%last point is a special case because there is no point after it, => no
%displacement and displacement vector
i = vecLength;

lengthI = mtLength(i);
simData(i).coord = [0 0 0; lengthI 0 0];
simData(i).centroid = [lengthI/2 0 0];
[simData(i).distanceMatrix, distanceVectorMatrix] = distMat(simData(i).coord); 
simData(i).distanceVectorMatrixN = ... 
    distanceVectorMatrix./repmat(simData(i).distanceMatrix,[1 1 3]);

simData(i).time = (i-1)*timeStep; 
simData(i).timepoint = i; 
simData(i).sigmaTime = aveInterval/2; 

simData(i).stats.qMatrix = diag([0.000064 0 0 mtLengthErr(i)^2 0 0 ...
        (0.000064+mtLengthErr(i)^2)/8 0 0]); 
simData(i).stats.noise = [1 1];

inputDataEntry.anaDat = simData;
inputDataEntry.fileInfo.identifier = 'GENERATED';
inputDataEntry.fileInfo.secondPath = 'by Monte Carlo simulations';

warning on MATLAB:divideByZero;
