function [arParam,maParam,noiseSigma,trajP,errFlag] = mlestim(traj,arOrder,maOrder,arParam0,maParam0)
%MLESTIM estimates the parameters of an ARMA modelto reproduce an observed trajectory
%
%SYNOPSIS [arParam,maParam,noiseSigma,trajP,errFlag] = mlestim(traj,arOrder,maOrder,arParam0,maParam0)
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
if arOrder < 0
    disp('--mlestim: "arOrder" should be a nonnegative integer!');
    errFlag = 1;
end
if maOrder < 0
    disp('--mlestim: "maOrder" should be a nonnegative integer!');
    errFlag = 1;
end
if arOrder == 0 && maOrder == 0
    disp('--mlestim: Either "arOrder" or "maOrder" should be nonzero!');
    errFlag = 1;
end
if errFlag
    disp('--mlestim: Please fix input data!');
    arParam = [];
    maParam = [];
    noiseSigma = [];
    return
end
if arOrder ~= 0
    [nRow,nCol] = size(arParam0);
    if nRow ~= 1
        disp('--innovPredict: "arParam0" should be a row vector!');
        errFlag = 1;
    else
        if nCol ~= arOrder
            disp('--innovPredict: Wrong length of "arParam0"!');
            errFlag = 1;
        end
        r = abs(roots([-arParam0(end:-1:1) 1]));
        if ~isempty(find(r<=1))
            disp('--innovPredict: Causality requires the polynomial defining the autoregressive part of the model not to have any zeros for z <= 1!');
            errFlag = 1;
        end
    end
end
if maOrder ~= 0
    [nRow,nCol] = size(maParam0);
    if nRow ~= 1
        disp('--innovPredict: "maParam0" should be a row vector!');
        errFlag = 1;
    else
        if nCol ~= maOrder
            disp('--innovPredict: Wrong length of "maParam0"!');
            errFlag = 1;
        end
        r = abs(roots([maParam0(end:-1:1) 1]));
        if ~isempty(find(r<=1))
            disp('--innovPredict: Invertibility requires the polynomial defining the moving average part of the model not to have any zeros for z <= 1!');
            errFlag = 1;
        end
    end
end
if errFlag
    disp('--mlestim: Please fix input data!');
    arParam = [];
    maParam = [];
    noiseSigma = [];
    return
end

%initial guess of parameters
param0 = [arParam0 maParam0];

%minimize the reduced likelihood (defined in Eq. 5.2.12 of "Introduction to Time
%Series and Forecasting" by Brockwell and Davis) to get best set of parameters.
[params,minFunc,exitFlag,output] = fmincon(@redLikelihood,param0,A,b,Aeq,beq,...
    lb,ub,nonlcon,options,OTHER PARAM);

%assign parameters obtained through minimization
arParam = params(1:arOrder);
maParam = params(arOrder+1:end);

%get maximum likelihood linear prediction of trajectory
[trajP,innovCoef,innovErr,errFlag] = innovPredict(traj,arOrder,maOrder,arParam,maParam);
if errFlag
    noiseSigma = [];
    return
end

%get standard deviation of white noise in process
trajDiff = (trajP-traj).^2./innovErr(1:end-1);
noiseSigma = sqrt(mean(trajDiff));
