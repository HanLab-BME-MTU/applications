function redLikeliV = redLikelihood(param,arOrder,maOrder,traj)
%REDLIKELIHOOD calculates the reduced likelihood of the fit of an ARMA model to observed data
%
%SYNOPSIS redLikeliV = redLikelihood(param,arOrder,maOrder,traj)
%
%INPUT  param   : Set of parameters in ARMA model (concat. of arParam and maParam)
%       arOrder : Order of autoregressive part of process.
%       maOrder : Order of moving average part of process.
%       traj    : Observed trajectory.    
%
%OUTPUT redLikeliV: Value of reduced likelihood.
%
%Khuloud Jaqaman, February 2004

%check if correct number of arguments were used when function was called
if nargin ~= nargin('redLikelihood')
    disp('--redLikelihood: Incorrect number of input arguments!');
    errFlag  = 1;
    xPredicted = [];
    innovCoef = [];
    innovErr = [];
    return
end

%distribute parameters
arParam = param(1:arOrder);
maParam = param(arOrder+1:end);

%get 1-step predictions of trajectory using innovations algorithm with
%the current parameters
[trajP,innovCoef,innovErr,errFlag] = innovPredict(traj,...
    arOrder,maOrder,arParam,maParam,0);
if errFlag
    error('redLikeliV: Could not predict trajectory!');
end

% %mean prediction error
% meanError = mean((trajP-traj).^2./innovErr(1:end-1));
% 
% %reduced likelihood
% redLikeliV = log(meanError) + mean(log(innovErr(1:end-1)));

%relative square prediction error per time point
relError = (trajP-traj).^2./innovErr(1:end-1);

%reduced likelihood
% if maOrder ~= 0
    redLikeliV = log(mean(relError)) + mean(log(innovErr(1:end-1)));
% else
%     errLength = length(relError)
%     redLikeliV = log(sum(relError(1:arOrder))/errLength) + ...
%         sum(log(innovErr(1:arOrder))/errLength;
% end
