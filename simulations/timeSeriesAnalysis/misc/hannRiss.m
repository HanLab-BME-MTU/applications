function [arParam,maParam,residuals,noiseSigma,errFlag] = hannRiss(traj,arOrder,maOrder)
%HANNRISS estimates the parameters of an ARMA process using the Hannan-Rissanen algorithm
%
%SYNOPSIS [arParam,maParam,residuals,noiseSigma,errFlag] = hannRiss(traj,arOrder,maOrder)
%
%INPUT  traj   : Trajectory to be modeled.
%       arOrder: Order of AR part of proposed model.
%       maOrder: Order of MA part of proposed model
%
%OUTPUT arParam   : Estimated AR parameters in model.
%       maParam   : Estimated MA parameters in model.
%       residuals : Estimation residuals.
%       noiseSigma: Estimated standard deviation of white noise.
%       errFlag   : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, February 2004

errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= nargin('hannRiss')
    disp('--hannRiss: Incorrect number of input arguments!');
    errFlag  = 1;
    arParam = [];
    return
end

%check input data
if arOrder < 0
    disp('--hannRiss: Variable "arOrder" should be >= 0!');
    errFlag = 1;
end
if maOrder < 1
    disp('--hannRiss: Variable "maOrder" should be >= 1!');
    errFlag = 1;
end
maxOrder = max(arOrder,maOrder);
trajLength = length(traj);
if trajLength < 25*maxOrder
    disp('--hannRiss: Length of trajectory should be at least 25 times larger than maximum model order!');
    errFlag = 1;
end
if errFlag
    disp('--hannRiss: please fix input data!');
    return
end

%STEP 1: Estimate noise by fitting to a high-order AR model

%assign order of initial AR model
arOrderI = 100;

%get parameters of initial AR model
[arParamI,noiseSigmaI,errFlag] = yuleWalker(traj,arOrderI);

%previous points to be used in AR prediction
prevPoints = zeros(trajLength-arOrderI,arOrderI);
for i=1:trajLength-arOrderI
    prevPoints(i,:) = traj(i+arOrderI-1:-1:i)';
end

%estimate noise at every time point
noise = zeros(trajLength,1);
noise(arOrderI+1:end) = traj(arOrderI+1:end) - prevPoints*arParamI';

%STEP 2: Estimate ARMA parameters by least squares linear regression

%previous points to be used in ARMA predicion
combOrder = arOrderI+maOrder;
prevPoints = zeros(trajLength-combOrder,arOrder+maOrder);
for i=1:trajLength-combOrder
    prevPoints(i,:) = [traj(i+combOrder-1:-1:i+combOrder-arOrder)' ...
            noise(i+combOrder-1:-1:i+arOrderI)'];
end

%get parameters (See p. 157 in Introduction to Time Series and Forecasting)
param = (prevPoints'*prevPoints)\(prevPoints'*traj(combOrder+1:end));

%distribute parameters
arParam = param(1:arOrder)';
maParam = param(arOrder+1:end)';

%get residuals
residuals = traj(combOrder+1:end) - prevPoints*param;

%compute sum of squared residuals
sqErr = sum(residuals.^2);

%estimate of white noise standard deviation
noiseSigma = sqrt(sqErr/(trajLength-combOrder));

%Omitted Step 3 (p.157), which is meant to improve the answer obtained in Step 2.
