function [arParamK,maParamK,arParamL,maParamL,varCovMat,wnVariance,...
    wnVector,aic,errFlag] = armaIdentKalman(trajectories,modelParam,minOpt)
%ARMAIDENTKALMAN determines the ARMA model (i.e. its order, coefficients and white noise variance) that best fits a time series which could have missing data points.
%
%SYNOPSIS [arParamK,maParamK,arParamL,maParamL,varCovMat,wnVariance,...
%    wnVector,aic,errFlag] = armaIdentKalman(trajectories,modelParam,minOpt)
%
%INPUT  trajectories: Observations of time series to be fitted. Either an 
%                     array of structures traj(1:nTraj).observations, or a
%                     2D array representing one single trajectory. 
%           .observations: 2D array of measurements and their uncertainties.
%                     Missing points should be indicated with NaN.
%       modelParam  : Structure array of models to be tested. Each entry consists of 
%                    2 elements:
%           .arParamP0   : Initial guess of partial autoregressive coefficients (row vector).
%           .maParamP0   : Initial guess of partial moving average coefficients (row vector).
%       minOpt      : Optional. Minimization option: 
%                     -'ml' for Matlab local minimizer "fmincon";
%                     -'tl' for Tomlab local minimizer "ucSolve";
%                     -'tg' for Tomlab global minimizer "glbFast"' followed
%                       by Tomlab local minimizer "ucSolve";
%                     -'nag' for NAG's local minimizerE04JAF.
%                     Default: 'ml'
%
%OUTPUT arParamK    : Estimated AR parameters using likelihood maximization.
%       maParamK    : Estimated MA parameters using likelihood maximization.
%       arParamL    : Estimated AR parameters using least squares fitting.
%       maParamL    : Estimated MA parameters using least squares fitting.
%       varCovMat   : Variance-covariance matrix of ARMA coefficients, 
%                     estimated via least squares fitting.
%       wnVariance  : Estimated variance of white noise in process.
%       wnVector    : Structure array containing the field:
%           .observations: Estimated white noise series in corresponding 
%                          trajectory.
%       aic         : Akaike's Information Criterion.
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%REMARKS Data is fitted using "armaCoefKalman" beginning with the different
%        input models. The results from the models are compared and the 
%        result that minimizes Akaike's Information Criterion is given as 
%        the final answer.
%
%Khuloud Jaqaman, July 2004

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arParamK = [];
maParamK = [];
arParamL = [];
maParamL = [];
varCovMat = [];
wnVariance = [];
wnVector = [];
aic = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
%all other validation is done in armaCoefKalman
if nargin < 2
    disp('--armaIdentKalman: Incorrect number of input arguments!');
    errFlag = 1;
    return
end

%check whether minOpt has one of the required values
if nargin == 2
    minOpt = 'ml';
else
    if (~strcmp(minOpt,'ml') && ~strcmp(minOpt,'tl') ...
            && ~strcmp(minOpt,'tg') && ~strcmp(minOpt,'nag'))
        disp('--armaIdentKalman: "minOpt" should be either "ml", "tl", "tg" or "nag"!');
        errFlag = 1;
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Comparison of ARMA models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%assign initial value of AIC
aic = 1e10; %(ridiculously large number)

%go over all suggested models
for i = 1:length(modelParam)
    
    %assign model data
    arParamP0 = modelParam(i).arParamP0;
    maParamP0 = modelParam(i).maParamP0;
    
    %estimate ARMA coeffients and white noise variance
    try
        [arParamK1,maParamK1,arParamL1,maParamL1,varCovMat1,wnVariance1,...
            wnVector1,aic1,errFlag] = armaCoefKalman(trajectories,...
            arParamP0,maParamP0,minOpt);
        if errFlag
            aic1 = [];
        end
    catch
        disp(['--armaIdentKalman: armaCoefkalman could not be executed for model #' num2str(i)]);
        aic1 = [];
    end
    
    %look for model with minimum AIC
    if isreal(aic1) %do not consider models that yield a complex AIC 
                    %(i.e. where minimization messed up)
        if aic1 < aic %if current AIC is smaller than minimum AIC so far
            arParamK = arParamK1; %update model
            maParamK = maParamK1;
            arParamL = arParamL1;
            maParamL = maParamL1;
            varCovMat = varCovMat1;
            wnVariance = wnVariance1;
            wnVector = wnVector1;
            aic = aic1;
        end
    else
        disp(['complex aic for model #' num2str(i) '!']);
    end %(if isreal(aic1))

end %(for i = 1:length(modelParam))


%%%%% ~~ the end ~~ %%%%%
