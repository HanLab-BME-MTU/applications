function [arParam,noiseSigma,arParamSigma,errFlag] = arlsestim0(traj,arOrder)
%ARLSESTIM0 estimates parameters of an AR model using least square fitting
%
%SYNOPSIS [arParam,noiseSigma,arParamSigma,errFlag] = arlsestim0(traj,arOrder)
%
%INPUT  traj   : Trajectory to be modeled (with measurement uncertainty.
%       arOrder: Order of proposed AR model.
%
%OUTPUT arParam     : Estimated parameters in model.
%       noiseSigma  : Estimated standard deviation of white noise.
%       arParamSIgma: Uncertainty in estimated AR coefficients.
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, February 2004

errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= nargin('arlsestim0')
    disp('--arlsestim0: Incorrect number of input arguments!');
    errFlag  = 1;
    arParam = [];
    noiseSigma = [];
    arParamSigma = [];
    return
end

%check input data
if arOrder < 1
    disp('--arlsestim0: Variable "arOrder" should be >= 1!');
    errFlag = 1;
end
[trajLength,nCol] = size(traj);
if trajLength < 5*arOrder
    disp('--arlsestim0: Length of trajectory should be at least 5 times larger than model order!');
    errFlag = 1;
end
if nCol ~= 2
    disp('--arlsestim0: "traj" should have one column for measurement and one for measurement uncertainty!');
    errFlag = 1;
end
if errFlag
    disp('--arlsestim0: please fix input data!');
    arParam = [];
    noiseSigma = [];
    arParamSigma = [];
    return
end
    
%previous points to be used in AR prediction 
%[(trajLength-arOrder) by arOrder matrix]
prevPoints = zeros(trajLength-arOrder,arOrder);
for i=1:trajLength-arOrder
    prevPoints(i,:) = traj(i+arOrder-1:-1:i,1)';
end

%weights matrix obtained from measurement uncertainties
%[(trajLength-arOrder) by (trajLength-arOrder) matrix]
weights = diag(traj(arOrder+1:end,2).^(-2));

%multiply weights by transpose of prevPoints 
%[arOrder by (trajLength-arOrder) matrix]
weightsP = prevPoints'*weights;

%get right hand side of set of equation 
%[arOrder by 1 vector]
rhs = weightsP*traj(arOrder+1:end,1);

%get matrix multiplying unknown parameters in left hand side of set of equation 
%[arOrder by arOrder matrix]
lhs = weightsP*prevPoints;

%clear variable
clear weightsP;

%variance-covariance matrix
varCovarMat = inv(lhs);

%solution
arParam = (varCovarMat*rhs)';

%check for causality of estimated model
r = abs(roots([-arParam(end:-1:1) 1]));
if ~isempty(find(r<=1.00001))
    disp('--arlsestim0: Warning: Predicted model not causal!');
end

%prediction error
epsilon = prevPoints*arParam' - traj(arOrder+1:end,1);

%standard deviation of white noise
noiseSigma = std(epsilon);

%variance-covariance matrix
arParamSigma = varCovarMat*sqrt((epsilon'*weights*epsilon)/(trajLength-2*arOrder));
