function [arParam,maParam,noiseSigma,trajP,aicc,errFlag] = mlestim(traj,...
    arOrder,maOrder,arParam0,maParam0)
%MLESTIM estimates the parameters of an ARMA model to reproduce an observed trajectory
%
%SYNOPSIS [arParam,maParam,noiseSigma,trajP,aicc,errFlag] = mlestim(traj,...
%    arOrder,maOrder,arParam0,maParam0)
%
%INPUT  traj    : Observed trajectory    
%       arOrder : Order of autoregressive part of process.
%       maOrder : Order of moving average part of process.
%       arParam0: Initial guess of AR coefficients.
%       maParam0: Initial guess of MA coefficients.
%
%OUTPUT arParam   : Coefficients of AR part.
%       maParam   : Coefficients of MA part.
%       noiseSigma: Standard deviation of white noise in model.
%       trajP     : Predicted trajectory using the innovations algorithm and
%                   the maximum likelihood estimates of the ARMA coefficients. 
%       aicc      : Bias-Corrected Akaike's Information Criterion for order selection.
%       errFlag   : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, February 2004

errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= nargin('mlestim')
    disp('--mlestim: Incorrect number of input arguments!');
    errFlag  = 1;
    arParam = [];
    maParam = [];
    noiseSigma = [];
    return
end

%check input data
[trajLength,nCol] = size(traj);
if nCol ~= 1
    disp('--mlestim: Variable "traj" should be a column vector!');
    errFlag = 1;
end
if arOrder < 1
    disp('--mlestim: "arOrder" should be >= 1!');
    errFlag = 1;
end
if maOrder < 1
    disp('--mlestim: "maOrder" should be >= 1!');
    errFlag = 1;
end
if errFlag
    disp('--mlestim: Please fix input data!');
    arParam = [];
    maParam = [];
    noiseSigma = [];
    return
end
[nRow,nCol] = size(arParam0);
if nRow ~= 1
    disp('--mlestim: "arParam0" should be a row vector!');
    errFlag = 1;
else
    if nCol ~= arOrder
        disp('--mlestim: Wrong length of "arParam0"!');
        errFlag = 1;
    end
    r = abs(roots([-arParam0(end:-1:1) 1]));
    if ~isempty(find(r<=1.00001))
        disp('--mlestim: Causality requires the polynomial defining the autoregressive part of the model not to have any zeros for z <= 1!');
        errFlag = 1;
    end
end
[nRow,nCol] = size(maParam0);
if nRow ~= 1
    disp('--mlestim: "maParam0" should be a row vector!');
    errFlag = 1;
else
    if nCol ~= maOrder
        disp('--mlestim: Wrong length of "maParam0"!');
        errFlag = 1;
    end
    r = abs(roots([maParam0(end:-1:1) 1]));
    if ~isempty(find(r<=1.00001))
        disp('--mlestim: Invertibility requires the polynomial defining the moving average part of the model not to have any zeros for z <= 1!');
        errFlag = 1;
    end
end
if errFlag
    disp('--mlestim: Please fix input data!');
    arParam = [];
    maParam = [];
    noiseSigma = [];
    return
end

%initial set of parameters
param0 = [arParam0 maParam0];
pLength = length(param0);

%define optimization options.
options = optimset('Display','iter');

%minimize the reduced likelihood (defined in Eq. 5.2.12 of "Introduction to Time
%Series and Forecasting" by Brockwell and Davis) to get best set of parameters.
[params,minFunc,exitFlag,output] = fmincon(@mlestimObj,param0,[],[],[],[],...
    -10*ones(pLength,1),10*ones(pLength,1),@mlestimConst,options,arOrder,maOrder,traj);

%assign parameters obtained through minimization
arParam = params(1:arOrder);
maParam = params(arOrder+1:end);

%check for causality and invertibility of estimated model
r = abs(roots([-arParam(end:-1:1) 1]));
if ~isempty(find(r<=1.00001))
    disp('--mlestim: Warning: Predicted model not causal!');
    errFlag = 1;
    noiseSigma = [];
    return
end
r = abs(roots([maParam(end:-1:1) 1]));
if ~isempty(find(r<=1.00001))
    disp('--mlestim: Warning: Predicted model not invertible!');
    errFlag = 1;
    noiseSigma = [];
    return
end
 
%If estimated model is OK, get maximum likelihood linear prediction of trajectory
[trajP,innovCoef,innovErr,errFlag] = innovPredict(traj,arOrder,maOrder,...
    arParam,maParam,1);
if errFlag
    noiseSigma = [];
    return
end

%get standard deviation of white noise in process
relError = (trajP-traj).^2./innovErr(1:end-1);
noiseSigma = sqrt(mean(relError));

%get Akaike's Information Criterion of model
aicc = trajLength*log(2*pi*noiseSigma^2) + sum(log(innovErr(1:end-1))) + ...
    trajLength + 2*trajLength*(arOrder+maOrder+1)/(trajLength-arOrder-maOrder-2);
