function [arParam,trajP,noiseSigma,errFlag] = arlsestim1(traj,arOrder,arParam0,noiseSigma0,multi)
%ARLSESTIM1 estimates parameters of an AR model when there are missing points
%
%SYNOPSIS [arParam,trajP,noiseSigma,errFlag] = arlsestim1(traj,arOrder,arParam0,noiseSigma0,multi)
%
%INPUT  traj        : Trajectory to be modeled (with measurement uncertainty).
%                     Missing points should be indicated with Inf.
%       arOrder     : Order of proposed AR model.
%       arParam0    : Initial value of parameters.
%       noiseSigma0 : Initial guess for white noise standard deviation.
%       multi       : 0/1 if fmincon/fminimax are to be used, respectively.
%
%OUTPUT arParam   : Estimated parameters in model.
%       trajP     : Measurement error-free predicted trajectory.
%       noiseSigma: Estimated standard deviation of white noise.
%       errFlag   : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, February 2004

errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= nargin('arlsestim1')
    disp('--arlsestim1: Incorrect number of input arguments!');
    errFlag  = 1;
    arParam = [];
    return
end

%check input data
if arOrder < 1
    disp('--arlsestim1: Variable "arOrder" should be >= 1!');
    errFlag = 1;
end
[trajLength,nCol] = size(traj);
if trajLength < 5*arOrder
    disp('--arlsestim1: Length of trajectory should be at least 5 times larger than model order!');
    errFlag = 1;
end
if nCol ~= 2
    disp('--arlsestim1: "traj" should have one column for measurement and one for measurement uncertainty!');
    errFlag = 1;
end
[nRow,nCol] = size(arParam0);
if nRow ~= 1
    disp('--arlsestim1: "arParam0" should be a row vector!');
    errFlag = 1;
else
    if nCol ~= arOrder
        disp('--arlsestim1: Wrong length of "arParam0"!');
        errFlag = 1;
    end
    r = abs(roots([-arParam0(end:-1:1) 1]));
    if ~isempty(find(r<=1.00001))
        disp('--arlsestim1: Causality requires the polynomial defining the autoregressive part of the model not to have any zeros for z <= 1!');
        errFlag = 1;
    end
end
if noiseSigma0 < 0
    disp('--arlsestim1: White noise standard deviation should be nonnegative!');
    errFlag = 1;
end
if multi~=0 && multi~=1
    disp('--arlsestim1: Variable "multi" should be 0 or 1!');
    errFlag = 1;
end
if errFlag
    disp('--arlsestim1: please fix input data!');
    return
end
    
%initial set of AR parameters and error-free measurements
unknown0 = [arParam0'; traj(:,1)];
indx = find(unknown0 == Inf); %find missing points
if ~isempty(indx)
    indxLow = indx(find(indx <= 2*arOrder)); %missing points at times <= arOrder
    indx = indx(find(indx > 2*arOrder)); %missing points at times > arOrder
    if indxLow(1) == arOrder+1
        unknown0(arOrder+1) = 0;
        indxLow = indxLow(2:end);
    end
    for i = indxLow %fill missing points at times <= arOrder
        unknown0(i) = arParam0(1:i-arOrder-1)*unknown0(i-1:-1:arOrder+1);
    end
    for i = indx %fill missing points at times > arOrder
        unknown0(i) = arParam0*unknown0(i-1:-1:i-arOrder);
    end
end

%find indices of available points
indx = find(traj(:,1)~=Inf);

%lower limit of measurement error-free observation
minVal = Inf*ones(trajLength,1);
minVal(indx) = traj(indx,1) - 3*traj(indx,2);
dummy = min(minVal([indx]));
minVal(find(minVal==Inf)) = dummy; %make the lower limit of missing point the minimum lower limit of observed points

%upper limit of measurement error-free observation
maxVal = -Inf*ones(trajLength,1);
maxVal(indx) = traj(indx,1) + 3*traj(indx,2);
dummy = max(maxVal([indx]));
maxVal(find(maxVal==-Inf)) = dummy; %make the upper limit of missing points the maximum upper limit of observed points

%define optimization options.
options = optimset('Display','iter','MaxFunEvals',100000,'maxIter',1000,...
    'TolCon',1e-8,'TolFun',1e-8,'TolX',1e-8);
%,'DerivativeCheck','on');

%minimize the sum of square errors to get best set of parameters.
if multi
    unknowns = fminimax(@arlsestim1Obj,unknown0,[],[],[],[],...
        [-10*ones(arOrder,1); minVal],[10*ones(arOrder,1); maxVal],...
        @arlsestim1Const,options,arOrder,traj,noiseSigma0,multi);
else
    unknowns = fmincon(@arlsestim1Obj,unknown0,[],[],[],[],...
        [-10*ones(arOrder,1); minVal],[10*ones(arOrder,1); maxVal],...
        @arlsestim1Const,options,arOrder,traj,noiseSigma0,multi);
end

%assign parameters obtained through minimization
arParam = unknowns(1:arOrder)';
trajP = unknowns(arOrder+1:end);

%check for causality of estimated model
r = abs(roots([-arParam(end:-1:1) 1]));
if ~isempty(find(r<=1.00001))
    disp('--arlsestim1: Warning: Predicted model not causal!');
end

%get standard deviation of white noise
prevPoints = zeros(trajLength-arOrder,arOrder); %previous points to be used in AR prediction
for i=1:trajLength-arOrder
    prevPoints(i,:) = trajP(i+arOrder-1:-1:i,1)';
end
noiseSigma = std(trajP(arOrder+1:end,1)-prevPoints*arParam'); %standard deviation
