function [xPredicted,innovCoef,innovErr,errFlag] = innovPredict(xObserved,...
    arOrder,maOrder,arParam,maParam,checkRoots)
%INNOVPREDICT calculates the 1-step predictor of a time series assuming it's an ARMA process
%
%SYNOPSIS [xPredicted,innovCoef,innovErr,errFlag] = innovPredict(xObserved,...
%    arOrder,maOrder,arParam,maParam,checkRoots)
%
%INPUT  xObserved   : Observed time series.
%       arOrder     : Order of autoregressive part of process.
%       maOrder     : Order of moving average part of process.
%       arParam     : Row vector of autoregression coefficients.
%       maPAram     : Row vector of moving average coefficients.
%       checkRoots  : 1 if AR and MA polynomial roots are to be checked, 0 otherwise.
%
%OUTPUT xPredicted  : Predicted time series, using the innovations
%                     algorithm and 1-step prediction.
%       innovCoef   : Coefficients calculated in algorithm.
%       innovErr    : Mean squared error of predictions.
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, February 2004

errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= nargin('innovPredict')
    disp('--innovPredict: Incorrect number of input arguments!');
    errFlag  = 1;
    xPredicted = [];
    innovCoef = [];
    innovErr = [];
    return
end

%check input data
if arOrder < 1
    disp('--innovPredict: "arOrder" should be >= 1!');
    errFlag = 1;
end
if maOrder < 1
    disp('--innovPredict: "maOrder" should be >= 1!');
    errFlag = 1;
end
if errFlag
    disp('--innovPredict: Please fix input data!');
    xPredicted = [];
    innovCoef = [];
    innovErr = [];
    return
end
[nRow,nCol] = size(arParam);
if nRow ~= 1
    disp('--innovPredict: "arParam" should be a row vector!');
    errFlag = 1;
else
    if nCol ~= arOrder
        disp('--innovPredict: Wrong length of "arParam"!');
        errFlag = 1;
    end
    if checkRoots
        r = abs(roots([-arParam(end:-1:1) 1]));
        if ~isempty(find(r<=1.00001))
            disp('--innovPredict: Causality requires the polynomial defining the autoregressive part of the model not to have any zeros for z <= 1!');
            errFlag = 1;
        end
    end
end
[nRow,nCol] = size(maParam);
if nRow ~= 1
    disp('--innovPredict: "maParam" should be a row vector!');
    errFlag = 1;
else
    if nCol ~= maOrder
        disp('--innovPredict: Wrong length of "maParam"!');
        errFlag = 1;
    end
    if checkRoots
        r = abs(roots([maParam(end:-1:1) 1]));
        if ~isempty(find(r<=1.00001))
            disp('--innovPredict: Invertibility requires the polynomial defining the moving average part of the model not to have any zeros for z <= 1!');
            errFlag = 1;
        end
    end
end
if errFlag
    disp('--innovPredict: Please fix input data!');
    xPredicted = [];
    innovCoef = [];
    innovErr = [];
    return
end

trajLength = length(xObserved); %trajectory length

%get autocovariance matrix of transformed process to use in innovations algorithm
[kappa,errFlag] = transAutoCov(arOrder,maOrder,arParam,maParam,...
    trajLength+2,trajLength+2,checkRoots);
if errFlag
    xPredicted = [];
    innovCoef = [];
    innovErr = [];
    return
end

%get innovations algorithm coefficients and prediction errors 
[innovErr,innovCoef,errFlag] = innovationsAlg(kappa,trajLength);
if errFlag
    xPredicted = [];
    innovCoef = [];
    innovErr = [];
    return
end

maxOrder = max(arOrder,maOrder); %needed for prediction

%initialize xPredicted
xPredicted = zeros(trajLength,1); %(not that, by definition, xPredicted(0) = 0)

%predict values between n=2 and n=maxOrder-1
for n=1:maxOrder-1
    xPredicted(n+1) = innovCoef(n,1:n)*(xObserved(n:-1:1)-xPredicted(n:-1:1));
end

%predict values between n=maxOrder and n=trajLength
for n=maxOrder:trajLength-1
    xPredicted(n+1) = arParam*xObserved(n:-1:n+1-arOrder) + innovCoef(n,1:maOrder)*...
        (xObserved(n:-1:n+1-maOrder)-xPredicted(n:-1:n+1-maOrder));
end
