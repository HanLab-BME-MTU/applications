function [arParam,noiseSigma,errFlag] = yuleWalker(traj,arOrder)
%YULEWALKER estimates the parameters of an AR process using the Yule-Walker Equations
%
%SYNOPSIS [arParam,noiseSigma,errFlag] = yuleWalker(traj,arOrder)
%
%INPUT  traj   : Trajectory to be modeled.
%       arOrder: Order of proposed AR model.
%
%OUTPUT arParam   : Estimated parameters in model.
%       noiseSigma: Estimated standard deviation of white noise.
%       errFlag   : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, February 2004

errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= nargin('yuleWalker')
    disp('--yuleWalker: Incorrect number of input arguments!');
    errFlag  = 1;
    arParam = [];
    return
end

%check input data
if arOrder < 1
    disp('--yuleWalker: Variable "arOrder" should be >= 1!');
    errFlag = 1;
end
if length(traj) < 5*arOrder
    disp('--yuleWalker: Length of trajectory should be at least 5 times larger than model order!');
    errFlag = 1;
end
if errFlag
    disp('--yuleWalker: please fix input data!');
    return
end
    
%get autocovariance function of trajectory
gamma = xcov(traj,arOrder,'unbiased');

%retain terms for nonnegative lags only
gamma = gamma(arOrder+1:end);

%form covariance matrix
covMat = diag(gamma(1)*ones(arOrder,1));
for i = 1:arOrder-1
    entry = gamma(i+1)*ones(arOrder-i,1);
    covMat = covMat + diag(entry,i) + diag(entry,-i);
end

%get parameters
arParam = (covMat\gamma(2:end))';

%get standard deviation of white noise
noiseSigma = sqrt(gamma(1)-arParam*gamma(2:end));
