function [arParam,noiseSigma,arParamSigma,errFlag] = arlsestim2(traj,arOrder,arParam0,noiseSigma0,arTol,sigmaTol)
%ARLSESTIM2 estimates parameters of an AR model using least square fitting, accounting for missing data points
%
%SYNOPSIS [arParam,noiseSigma,arParamSigma,errFlag] = arlsestim2(traj,arOrder,arParam0)
%
%INPUT  traj        : Trajectory to be modeled (with measurement uncertainty).
%                     Missing points should be indicated with NaN.
%       arOrder     : Order of proposed AR model.
%       arParam0    : Initial guess of AR coefficients.
%       noiseSigma0 : Initial guess of white noise standard deviation.
%       arTol       : Change in AR coefficients that below which it is safe to stop.
%       sigmaTol    : Change in WN standard deviation below which it is safe to stop.
%
%OUTPUT arParam     : Estimated parameters in model.
%       noiseSigma  : Estimated standard deviation of white noise.
%       arParamSigma: Uncertainty in estimated AR coefficients.
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, February 2004

errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= nargin('arlsestim2')
    disp('--arlsestim2: Incorrect number of input arguments!');
    errFlag  = 1;
    arParam = [];
    noiseSigma = [];
    return
end

%check input data
if arOrder < 1
    disp('--arlsestim2: Variable "arOrder" should be >= 1!');
    errFlag = 1;
end
[trajLength,nCol] = size(traj);
if trajLength < 5*arOrder
    disp('--arlsestim2: Length of trajectory should be at least 5 times larger than model order!');
    errFlag = 1;
end
if nCol ~= 2
    disp('--arlsestim2: "traj" should have one column for measurement and one for measurement uncertainty!');
    errFlag = 1;
end
[nRow,nCol] = size(arParam0);
if nRow ~= 1
    disp('--arlsestim2: "arParam0" should be a row vector!');
    errFlag = 1;
else
    if nCol ~= arOrder
        disp('--arlsestim2: Wrong length of "arParam0"!');
        errFlag = 1;
    end
    r = abs(roots([-arParam0(end:-1:1) 1]));
    if ~isempty(find(r<=1.00001))
        disp('--arlsestim2: Causality requires the polynomial defining the autoregressive part of the model not to have any zeros for z <= 1!');
        errFlag = 1;
    end
end
if noiseSigma0 < 0
    disp('--arlsestim2: White noise standard deviation should be nonnegative!');
    errFlag = 1;
end
if arTol <= 0
    disp('--arlsestim2: Variable "arTol" should be positive!');
    errFlag = 1;
end
if sigmaTol <= 0
    disp('--arlsestim2: Variable "sigmaTol" should be positive!');
    errFlag = 1;
end
if errFlag
    disp('--arlsestim2: please fix input data!');
    arParam = [];
    noiseSigma = [];
    return
end

indx = find(isnan(traj(:,1))); %find missing points

%check if there are missing points for time < arOrder. Algorithm cannot
%proceed if there are such points.
if ~isempty(indx)
    if indx(1) <= arOrder
        disp('--arlsestim2: There are missing points for time points < arOrder!');
        disp('  Please fix input data such that at least the first arOrder time points are available!');
        errFlag = 1;
        arParam = [];
        noiseSigma = [];
        return
    end
end

%"measurement uncertainty" that is assigned to missing points
maxErr = 1000*max(traj(:,2));
traj(indx,2) = maxErr;

%initialize unknowns
arParam = arParam0;
noiseSigma = noiseSigma0;

%store "traj" in new variable "trajNew"
trajNew = traj;

%changes in AR parameters and WNSD, to be compared with "arTol" and "sigmaTol" 
%to determine whether solution has been reached.
arParamDiff = 10*arTol;
sigmaDiff = 10*sigmaTol;

while arParamDiff > arTol || sigmaDiff > sigmaTol
    
    %predict missing points using proposed AR model
    for i = indx 
        trajNew(i,1) = arParam*trajNew(i-1:-1:i-arOrder,1);
    end
    
    %modify measurement uncertainty to include white noise standard deviation as well
    variance = noiseSigma^2;
    trajNew(:,2) = sqrt(traj(:,2).^2+variance);
    
    %call arlsestim0 to update/predict parameters and white noise standard deviation
    [arParamNew,noiseSigmaNew,arParamSigma,errFlag] = arlsestim0(trajNew,arOrder);
    
    %compute variables to be compared with "arTol" and "sigmaTol"
    arParamDiff = max(abs(arParamNew-arParam)); %maximum change in AR coefficients
    sigmaDiff = abs(noiseSigmaNew-noiseSigma); %change in white noise standard deviation
    
    %update unknowns
    arParam = arParamNew;
    noiseSigma = noiseSigmaNew;
    
end

%check for causality of estimated model
r = abs(roots([-arParam(end:-1:1) 1]));
if ~isempty(find(r<=1.00001))
    disp('--arlsestim2: Warning: Predicted model not causal!');
end

