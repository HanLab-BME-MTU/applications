function [sumSquareErr] = arlsestim1Obj(unknown,arOrder,traj,noiseSigma,multi)
%ARLSESTIM1OBJ computes the sum of square prediction errors for arlsestim1
%
%SYNOPSIS [sumSquareErr,sseGrad] = arlsestim1Obj(unknown0,arOrder,traj,noiseSigma,multi)
%
%INPUT  unknown     : AR parameters and measurement error-free trajectory.
%       arOrder     : Order of proposed AR model.
%       traj        : Observed trajectory to be modeled (with measurement uncertainty).
%                     Missing points should be indicated with Inf.
%       noiseSigma  : Standard deviation of white noise.
%       multi       : 0/1 if fmincon/fminimax are being used, respectively.
%
%OUTPUT sumSquareErr: sum of square prediction errors.
%       sseGrad     : vector of partial derivatives of sumSquareErr.
%
%Khuloud Jaqaman, February 2004

errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= nargin('arlsestim1Obj')
    disp('--arlsestim1Obj: Incorrect number of input arguments!');
    errFlag  = 1;
    arParam = [];
    return
end

%length of trajectory
trajLength = length(traj(:,1));

%check input data
if arOrder < 1
    disp('--arlsestim1Obj: Variable "arOrder" should be >= 1!');
    errFlag = 1;
end
[nRow,nCol] = size(unknown);
if nRow ~= trajLength+arOrder
    disp('--arlsestim1Obj: Variable "unknown" has the wrong length!');
    errFlag = 1;
end
if nCol ~= 1
    disp('--arlsestim1Obj: Variable "unknown" should be a column vector!');
    errFlag = 1;
end
if noiseSigma < 0
    disp('--arlsestim1Obj: White noise standard deviation should be nonnegative!');
    errFlag = 1;
end
if multi~=0 && multi~=1
    disp('--arlsestim1Obj: Variable "multi" should be 0 or 1!');
    errFlag = 1;
end
if errFlag
    disp('--arlsestim1Obj: please fix input data!');
    return
end

%AR coefficients
arParam = unknown(1:arOrder)';

%measurement error-free trajectory
trajNoErr = unknown(arOrder+1:end);

%find points where data is available
indx = find(traj(:,1) ~= Inf);
indxLow = indx(find(indx <= arOrder)); %points at times <= arOrder
indx = indx(find(indx > arOrder)); %points at times > arOrder

%trajectory length for times smaller than arOrder excluding missing points
indxLowLength = length(indxLow);

%trajectory length for times greater than arOrder excluding missing points
indxLength = length(indx);

%check for causality of model
r = abs(roots([-arParam(end:-1:1) 1]));

if ~isempty(find(r<=1.00001)) %if not causal 

    if multi %multi-objective
        sumSquareErr = 1e10*ones(2*indxLength+indxLowLength);
    else %single objective
        sumSquareErr = 1e10;
    end
%     sseGrad = [];
    
else %causal

    
    %map from original sequence for time>arOrder (with missing points) to sequence without missing points
    indxInv = zeros(trajLength-arOrder,1);
    indxInv(indx) = [1:indxLength];
    
    %set of observations
    trajWithErr = [traj(indx,1); traj(indxLow,1); traj(indx,1)];
    
    %previous points to be used in AR prediction xHat(t) + e(t) + mu(t) = sum_(i=1)^p[a_ix(t-i)]
    prevPoints = zeros(indxLength,arOrder);
    for i=1:indxLength
        j = indx(i);
        prevPoints(i,:) = trajNoErr(j-1:-1:j-arOrder)';
    end
    
    %get prediction errors
    errVec = zeros(2*indxLength+indxLowLength,1);
    errVec(1:indxLength) = prevPoints*arParam' - traj(indx,1); %from AR model
    errVec(indxLength+1:indxLength+indxLowLength) = trajNoErr(indxLow) - traj(indxLow); %from measurement-error free points
    errVec(indxLength+indxLowLength+1:end) = trajNoErr(indx) - traj(indx); %from measurement-error free points
    
    %weights obtained from measurement uncertainties and white noise variance
    weights = zeros(2*indxLength+indxLowLength,1);
    weights(1:indxLength) = 1./(noiseSigma^2+traj(indx,2).^2);
    weights(indxLength+1:indxLength+indxLowLength) = traj(indxLow,2).^(-2);
    weights(indxLength+indxLowLength+1:end) = traj(indx,2).^(-2);
    
    %value of function to be minimized 
    if multi %vector of weighted square errors (multi-objective)
        sumSquareErr = errVec.*weights.*errVec/2;
    else %sum of weighted square errors (single objective)
        sumSquareErr = errVec'*(weights.*errVec)/2;
    end
    
%     %initialize partial derivatives vector
%     sseGrad = zeros(size(unknown));
%     
%     %partial derivatives w.r.t. AR coefficients
%     dummy = weights(1:indxLength).*errVec(1:indxLength);
%     sseGrad(1:arOrder) = prevPoints'*dummy;
%     
%     %partial derivatives w.r.t. measurement error-free data points 
%     %AR part
%     for i=arOrder+1:length(sseGrad)-1
%         j = i - arOrder;
%         subIndx = indxInv(max(j+1,arOrder+1):min(j+arOrder,trajLength));
%         subIndx = subIndx(find(subIndx));
%         backShift = indx(subIndx) - j;
%         sseGrad(i) = arParam(backShift)*(weights(subIndx).*errVec(subIndx));
%     end
%     %additional part
%     sseGrad(arOrder+indxLow) = sseGrad(indxLow) + weights(indxLength+1:indxLength+indxLowLength)...
%         .*errVec(indxLength+1:indxLength+indxLowLength);
%     sseGrad(arOrder+indx) = sseGrad(indx) + weights(indxLength+indxLowLength+1:end)...
%         .*errVec(indxLength+indxLowLength+1:end);
    
end
