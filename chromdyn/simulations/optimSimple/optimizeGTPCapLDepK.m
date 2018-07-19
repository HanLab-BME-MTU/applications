%SCRIPT "optimizeGTPCapLDepK" calls a series of functions to calibrate the parameters in mtGTPCapLDepK
%
%Khuloud Jaqaman 11/03

%initial state of MT
initialState = struct('mtLength0',0.44,'capSize0',3,'unitConc',20);

%initialize parameters to be callibrated.
%k*[T] is used instead of k in some cases to make all rates in units of
%s^-1. This should be better for optimization.
numParam = 10;      %number of parameters to be optimized
param0(1) = 0.4*initialState.unitConc; %kTOnTFree*[T] (s^-1)
param0(2) = 0.2*initialState.unitConc; %addAmpT*[T] (s^-1)
param0(3) = 10;     %addWidT (micrometer^-1)
param0(4) = 0.85;   %addLenT (micrometer)
param0(5) = 0.2*initialState.unitConc; %kTOnD*[T] (s^-1)
param0(6) = 14;     %kDOffFree (s^-1)
param0(7) = 7;      %addAmpD (s^-1)
param0(8) = 10;     %addWidD (micrometer^-1)
param0(9) = 0.85;   %addLenD(micrometer)
param0(10) = 10;    %kHydrolysis (s^-1)
param0 = param0';

%experimental data used for calibration
numData  = 6;                %number of different characteristics used for optimization
expData(1,:) = [1.86 0.03];  %growth speed (micrometers/min)
expData(2,:) = [1.69 0.03];  %shrinkage speed (micrometers/min)
expData(3,:) = [0.27 0.1];   %catastrophe frequency (s^-1)
expData(4,:) = [0.27 0.1];   %rescue frequency (s^-1)
expData(5,:) = [0.48 0.01];  %minumum length (micrometers)
expData(6,:) = [1.12 0.02];  %maximum length (micrometers)

%define additional parameters to be passed on to function to be minimized
lMin = 0.2;         %minimum length of microtubule (micrometers)
lMax = 1.5;         %maximum length of microtubule (micrometers)

%run information. Note that simTimeStep will be modified based on the
%values of the rates. The rest are fixed.
runInfo = struct('maxNumSim',1,'totalTime',10000,'simTimeStep',0.3,...
    'timeEps',0.5,'expTimeStep',1,'aveInterval',0.6,'analyzeOpt',2);

%do not save trajectories or statistics
saveTraj = struct('saveOrNot',0,'fileName',[]);
saveStats = struct('saveOrNot',0,'fileName',[]);

%define lower and upper bounds of parameters
lb = [0 0 6/(lMax-lMin) lMin 0 0 0 6/(lMax-lMin) lMin 0]'; %lower
ub = [20 20 100 lMax 10 30 30 100 lMax 25]';         %upper

%write the linear constraints on the parameters in the form Ax=b, where x
%is the set of parameters.
%these constraints are
%       addAmpT - kTOnTFree <= 0
%       addAmpD - kDOffFree <= 0
A = [-1 1 0 0 0 0 0 0 0 0; 0 0 0 0 0 -1 1 0 0 0];
b = [0; 0];
%Note that nonlinear constraints are written in the function "nonlinConLDepK"

%%%complete run

%define optimization options.
options = optimset('DiffMinChange',0.1,'DiffMaxChange',0.2,'GradConstr','on',...
    'LargeScale','off','TolCon',0.0001,'TolFun',0.000001,'TolPCG',1,'TolX',1,...
    'Display','iter');

%now minimize objective function "objectiveGTPCapLDepK".
[param1,objectiveF1,exitFlag1] = fmincon(@objectiveGTPCapLDepK,param0,A,b,[],[],...
    lb,ub,@nonlinConLDepK,options,numParam,lMin,lMax,initialState,runInfo,...
    saveTraj,saveStats,numData,expData);

