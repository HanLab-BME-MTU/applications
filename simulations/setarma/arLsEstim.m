function [arParam,noiseSigma,errFlag] = arLsEstim(traj,arOrder)
%ARLSESTIM estimates parameters of an AR model using least square fitting
%
%SYNOPSIS [arParam,noiseSigma,errFlag] = arLsEstim(traj,arOrder)
%
%INPUT  traj   : Trajectory to be modeled (with measurement uncertainty.
%       arOrder: Order of proposed AR model.
%
%OUTPUT arParam   : Estimated parameters in model.
%       noiseSigma: Estimated standard deviation of white noise.
%       errFlag   : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, February 2004

errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= nargin('arLsEstim')
    disp('--arLsEstim: Incorrect number of input arguments!');
    errFlag  = 1;
    arParam = [];
    return
end

%check input data
if arOrder < 1
    disp('--arLsEstim: Variable "arOrder" should be >= 1!');
    errFlag = 1;
end
[trajLength,nCol] = size(traj);
if trajLength < 5*arOrder
    disp('--arLsEstim: Length of trajectory should be at least 5 times larger than model order!');
    errFlag = 1;
end
if nCol ~= 2
    disp('--arLsEstim: "traj" should have one column for measurement and one for measurement uncertainty!');
    errFlag = 1;
end
if errFlag
    disp('--arLsEstim: please fix input data!');
    return
end
    
%previous points to be used in AR prediction
prevPoints = zeros(trajLength-arOrder,arOrder);
for i=1:trajLength-arOrder
    prevPoints(i,:) = traj(i+arOrder-1:-1:i,1)';
end

%weights matrix obtained from measurement uncertainties
weights = diag(traj(arOrder+1:end,2).^(-2));

%multiply weights by transpose of prevPoints
weights = prevPoints'*weights;

%get solution
arParam = inv(weights*prevPoints)*weights*traj(arOrder+1:end,1);
arParam = arParam';

noiseSigma = std(traj(arOrder+1:end,1)-prevPoints*arParam');
