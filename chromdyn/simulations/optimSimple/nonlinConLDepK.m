function [c,ceq,GC,GCeq] = nonlinConLDepK(param,numParam,lMin,lMax,...
    initialState,runInfo,saveTraj,saveStats,numData,expData)
%NONLINCONCLDEPK evaluates the nonlinear contraints on parameters in mtGTPCapLDepK

%SYNOPSIS [c,ceq,GC,GCeq] = nonlinConLDepK(param,numParam,lMin,lMax,...
%    initialState,runInfo,saveTraj,saveStats,numData,expData)
%
%INPUT 
%
%OUTPUT
%
%COMMENTS: Many of these input parameters are not needed, but it seems they must be
%passed because of the structure of fmincon.

%check if correct number of arguments were used when function was called
if nargin ~= nargin('nonlinConLDepK')
    disp('--analyzeMtTrajectory: Incorrect number of input arguments!');
    errFlag  = 1;
    return;
end

errFlag = 0;

if lMin < 0
    disp('--nonlinConLDepK: lMin cannot be negative!');
    errFlag = 1;
end
if lMax <= lMin
    disp('--nonlinConLDepK: lMax must be greater than lMin!');
    errFlag = 1;
end
if errFlag
    error('--nonlinConLDepK: Please fix input data!');
end


c = [lMin + 3/param(3) - param(4);...
        -lMax + 3/param(3) + param(4);...
        lMin + 3/param(8) - param(9);...
        -lMax + 3/param(8) + param(9)];

ceq = [];

GC = [0 0 0 0;...
        0 0 0 0;...
        -3/(param(3)^2) -3/(param(3)^2) 0 0;...
        -1 1 0 0;...
        0 0 0 0;...
        0 0 0 0;...
        0 0 0 0;...
        0 0 -3/(param(8)^2) -3/(param(8)^2);...
        0 0 -1 1;...
        0 0 0 0];

GCeq = [];
