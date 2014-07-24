function [inputDataEntry,errFlag] = getAnaDat3(mtLength,mtLengthErr,...
    simNum,timeStep,aveInterval)
%GETANADAT3 writes SPB to GFP Tag distance data in format required for calcMTDynamics. 

%SYNOPSIS [inputDataEntry,errFlag] = getAnaDat3(mtLength,mtLengthErr,...
%   simNum,timeStep,aveInterval)
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
%Khuloud Jaqaman 11/03

errFlag = 0;

vecLength = length(mtLength);

simData(1:vecLength)=struct('info',[],'coord',[],'centroid',[],...
    'distanceMatrix',[],'distanceVectorMatrixN',[],'displacementVectorN',[],...
    'displacement',[],'time',[],'timepoint',[],'sigmaTime',[],'stats',[]);

%general information
simData(1).info.name = sprintf('simMtDynamicsRun%i',simNum);
simData(1).info.idlisttype = 'synthetic';
simData(1).info.nTags = 2;
simData(1).info.labelColor = [{'spb'};{'tagPos'}];
simData(1).info.created = nowString;
simData(1).info.dataProperties = [];

warning off MATLAB:divideByZero;

%data
for i=1:vecLength-1
    
    lengthIx = mtLength(i,1); %x coordinate of tag at i (micrometers)
    lengthIy = mtLength(i,2); %y coordinate
    lengthIz = mtLength(i,3); %z coordinate
    xDisp = mtLength(i+1,1)-lengthIx; %displacement in each direction
    yDisp = mtLength(i+1,2)-lengthIy;
    zDisp = mtLength(i+1,3)-lengthIz;
    totalDisp = sqrt(xDisp^2+yDisp^2+zDisp^2); %magnitude of displacement
    if totalDisp ~= 0
        xDisp = xDisp/totalDisp; %unit vector in direction of displacement
        yDisp = yDisp/totalDisp;
        zDisp = zDisp/totalDisp;
    end
    
    simData(i).coord = [0 0 0; lengthIx lengthIy lengthIz]; %position of SPB and tag at i 
    %(micrometers)
    simData(i).centroid = [lengthIx/2 lengthIy/2 lengthIz/2]; %position of center of geometry (micrometers)
    [simData(i).distanceMatrix, distanceVectorMatrix] = ... %matrix of distances between  
        distMat(simData(i).coord); %SPB&SPB, SPB&MT-tip, MT-tip&SPB, MT-tip&MT-tip, and
    %matrices of vectors pointing from one to the other
    simData(i).distanceVectorMatrixN = ... %normalized "distanceVectorMatrix"
        distanceVectorMatrix./repmat(simData(i).distanceMatrix,[1 1 3]);
    
    simData(i).displacementVectorN = [0 0 0; xDisp yDisp zDisp];  %unit vector in direction of displacement
    simData(i).displacement = [0;totalDisp]; %magnitude of displacement between i and i+1 (micrometers)
    
    simData(i).time = (i-1)*timeStep; %time at iteration i (seconds)
    simData(i).timepoint = i; %iteration number, i.e.
    simData(i).sigmaTime = aveInterval/2; %interval of averaging (seconds) over 2
    
    simData(i).stats.qMatrix = diag([0.000064 0.000064 0.000064 ...   %variance in position of SPB, tag
            mtLengthErr(i,1)^2 mtLengthErr(i,2)^2 mtLengthErr(i,3)^2 ...  %and centroid (micrometers^2)
            (0.000064+mtLengthErr(i,1)^2)/8 (0.000064+mtLengthErr(i,2)^2)/8 ...
            (0.000064+mtLengthErr(i,3)^2)/8]); 
    simData(i).stats.noise = [1 1];
    
end

%last point is a special case because there is no point after it, => no
%displacement and displacement vector
i = vecLength;

lengthIx = mtLength(i,1);
lengthIy = mtLength(i,2);
lengthIz = mtLength(i,3);
simData(i).coord = [0 0 0; lengthIx lengthIy lengthIz];
simData(i).centroid = [lengthIx/2 lengthIy/2 lengthIz/2];
[simData(i).distanceMatrix, distanceVectorMatrix] = distMat(simData(i).coord); 
simData(i).distanceVectorMatrixN = ... 
    distanceVectorMatrix./repmat(simData(i).distanceMatrix,[1 1 3]);

simData(i).time = (i-1)*timeStep; 
simData(i).timepoint = i; 
simData(i).sigmaTime = aveInterval/2; 

simData(i).stats.qMatrix = diag([0.000064 0.000064 0.000064 ...
        mtLengthErr(i,1)^2 mtLengthErr(i,2)^2 mtLengthErr(i,3)^2 ...
        (0.000064+mtLengthErr(i,1)^2)/8 (0.000064+mtLengthErr(i,2)^2)/8 ...
        (0.000064+mtLengthErr(i,3)^2)/8]); 
simData(i).stats.noise = [1 1];

inputDataEntry.anaDat = simData;
inputDataEntry.fileInfo.identifier = 'GENERATED';
inputDataEntry.fileInfo.secondPath = 'by Monte Carlo simulations';

warning on MATLAB:divideByZero;
