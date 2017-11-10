function objective = objectiveGTPCapLDepK(param,numParam,lMin,lMax,...
    initialState,runInfo,saveTraj,saveStats,numData,expData)
%OBJECTIVEGTPCAPLDEPK evaluates the sum of squares of differences between characteristics of experimental trajectories and those simulated using mtGTPCapLDepK 
%
%SYNOPSIS objective = objectiveGTPCapLDepK(param,numParam,lMin,lMax,...
%   initialState,runInfo,numData,expData)
%
%INPUT  param : array of parameters to be optimized.
%       numParam : number of parameters to be optimized.
%       lMin  : minimum length of microtubule.
%       lMax  : maximum length of microtubule.
%       initialState : initial state of MT. See "analyzeMtTrajectory" for details.
%       runInfo  : simulation information. See "analyzeMtTrajectory" for details.
%       saveTraj : structure determining whether generated trajectories
%                  will be saved.
%       saveStats: structure determining whether statistical
%                  characteristics will be saved.
%       numData  : number of characterstics used in comparing.
%       expData  : values of characteristics from experiments.
%
%OUTPUT objective: value of the sum of squares of differences
%
%Khuloud Jaqaman, 11/03

%check if correct number of arguments were used when function was called
if nargin ~= nargin('objectiveGTPCapLDepK')
    error('--objectiveGTPCapLDepK: Incorrect number of input arguments!');
end

%check input data
errFlag = 0;
% if numParam ~= 10
%     disp('--objectiveGTPCapLDepK: "numParam" must equal 10!');
%     errFlag = 1;
% end
if numData ~= 6
    disp('--objectiveGTPCapLDepK: "numData" must equal 6!');
    errFlag = 1;
end
if errFlag
    error('--objectiveGTPCapLDepK: "numParam" and/or "numData" must be fixed!');
end
if lMin < 0
    disp('--objectiveGTPCapLDepK: lMin cannot be negative!');
    errFlag = 1;
end
if lMax <= lMin
    disp('--objectiveGTPCapLDepK: lMax must be greater than lMin!');
    errFlag = 1;
end
if length(param) ~= numParam
    disp('--objectiveGTPCapLDepK: Length of array "param" incorrect!');
    errFlag = 1;
end
if length(expData) ~= numData
    disp('--objectiveGTPCapLDepK: Length of array "expData" incorrect!');
    errFlag = 1;
end
if errFlag
    error('--objectiveGTPCapLDepK: Please fix input data!');
end

%define parameters to call "analyzeMtTrajectory"
concT = initialState.unitConc;
modelParam = struct('minLength',lMin,'maxLength',lMax,'kTOnTFree',param(1)/concT,...
    'addAmpT',param(2)/concT,'addWidT',param(3),'addLenT',param(4),'kTOff',0.0,...
    'kTOnD',param(5)/concT,'kDOffFree',param(6),'addAmpD',param(7),'addWidD',...
    param(8),'addLenD',param(9),'kHydrolysis',param(10));

kTOnTFreeEff = param(1); %effective rate constants of "unit" addition
kTOnDEff = param(5);

runInfo.simTimeStep = min(0.99*runInfo.timeEps/max([kTOnTFreeEff kTOnDEff...
        param(6) param(10)]),runInfo.aveInterval/2); %simulation time step 

%call analyzeMtTrajectory
[errFlag,dataStats] = analyzeMtTrajectory(2,modelParam,initialState,runInfo,...
    saveTraj,saveStats);
if errFlag
    error('--objectiveGTPCapLDepK: analyzeMtTrajectory did not finish successfully!');
end

%extract characterstics to be compared with experimental data
gSpeed = dataStats.growthSpeed;
if isempty(gSpeed)
    gSpeed(1) = Inf;
end
sSpeed = abs(dataStats.shrinkageSpeed);
if isempty(sSpeed)
    sSpeed(1) = Inf;
end
catFreq = dataStats.catFreq;
if isempty(catFreq)
    catFreq = Inf;
end
resFreq = dataStats.resFreq;
if isempty(resFreq)
    resFreq = Inf;
end
minLength = dataStats.minDistanceM5;
maxLength = dataStats.maxDistanceM5;

%calculate objective as the sum of the squares of the differences between
%characteristics of experimental and simulation trajectories.
objective = ((gSpeed(1)-expData(1,1))/expData(1,1))^2 + ((sSpeed(1)-expData(2,1))/expData(2,1))^2 ...
    + ((catFreq(1)-expData(3,1))/expData(3,1))^2 + ((resFreq(1)-expData(4,1))/expData(4,1))^2 ...
    + ((minLength(1)-expData(5,1))/expData(5,1))^2 + ((maxLength(1)-expData(6,1))/expData(6,1))^2;
%disp(sprintf('%f',objective));