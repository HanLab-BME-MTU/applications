function [arParamK,maParamK,wnVariance,wnVector,aic,varCovMat,arParamL,maParamL,errFlag] ...
    = armaIdentKalman(trajectories,modelParam)
%ARMAIDENTKALMAN uses Kalman prediction and filtering to determine the ARMA model (i.e. its order, coefficients and white noise variance) that best fits a time series which could have missing data points.
%
%SYNOPSIS [arParamK,maParamK,wnVariance,wnVector,aic,varCovMat,arParamL,maParamL,errFlag] ...
%    = armaIdentKalman(trajectories,modelParam)
%
%INPUT  trajectories: Structure array containing trajectories to be modeled:
%           .observations: 2D array of measurements and their uncertainties.
%                          Missing points should be indicated with NaN.
%       modelParam  : Structure array of models to be tested. Each entry consists of 
%                    2 elements:
%           .arParamP0   : Initial guess of partial autoregressive coefficients (row vector).
%           .maParamP0   : Initial guess of partial moving average coefficients (row vector).
%
%OUTPUT arParamK   : 1st row: AR coefficients estimated by Kalman filter.
%                    2nd row: corresponding partial AR coefficients.
%       maParamK   : 1st row: MA coefficients estimated by Kalman filter.
%                    2nd row: corresponding partial MA coefficients.
%       wnVariance : Estimated white noise variance.
%       wnVector   : Structure containing 1 field:
%           .series   : Estimated white noise series in a trajectory.
%       aic        : Akaike's Information Criterion.
%       varCovMat  : Variance-Covariance matrix of estimated ARMA coefficients.
%       arParamL   : AR parameters estimated by least square fitting.
%       maParamL   : MA parameters estimated by least square fitting.
%       errFlag    : 0 if function executes normally, 1 otherwise.
%
%REMARKS The ARMA model that best fits a group of time series is determined
%        using the Kalman filter techniques presented in R. H. Jones, "Maximum 
%        Likelihood Fitting of ARMA Models to Time Series with Missing Observations",
%        Technometrics 22: 389-395 (1980). This yields the ARMA
%        coefficients, white noise variance as well as white noise series
%        in the process. Using the latter, the problem is re-formulated as
%        a least square fitting problem, yielding an estimate of the
%        variance-covariance matrix of the predicted ARMA coefficients, as
%        well as another estimate of the ARMA coefficients. The latter
%        should agree with coefficients predicted using the Kalman filter!
%
%Khuloud Jaqaman, July 2004

%initialize output
arParamK = [];
maParamK = [];
wnVariance = [];
wnVector = [];
aic = [];
varCovMat = [];
arParamL = [];
maParamL = [];
errFlag = 0;

%check if correct number of arguments was used when function was called
if nargin < nargin('armaIdentKalman')
    disp('--armaIdentKalman: Incorrect number of input arguments!');
    errFlag = 1;
    return
end

%assign initial value of AIC
aic = 1e10; %(ridiculously large number)

%go over all suggested models
for i = 1:length(modelParam)
    
    %get model characteristics
    arParamP0 = modelParam(i).arParamP0;
    maParamP0 = modelParam(i).maParamP0;
    
    %estimate ARMA coeffients and white noise variance
    try
        [arParam1,maParam1,wnVariance1,wnVector1,aic1,aicc1,errFlag] = ...
            armaCoefKalman(trajectories,arParamP0,maParamP0);
        if errFlag
            aic1 = [];
            aicc1 = [];
        end
    catch
        disp('--armaIdentKalman: armaCoefkalman could not be executed for model #');
        disp(num2str(i));
        aic1 = [];
    end
    
    %compare current AIC to minimum AIC
    %if it is smaller, update results
    if isreal(aic1)
        if aic1 < aic
            arParamK = arParam1;
            maParamK = maParam1;
            wnVariance = wnVariance1;
            wnVector = wnVector1;
            aic = aic1;
        end
    else
        disp(['complex aic for model #' num2str(i) '!']);
    end

end %(for i = 1:length(modelParam))

%evaluate the variance-covariance matrix of the predicted ARMA coefficients
%by reformulating the problem as a least square fitting.
[varCovMat,arParamL,maParamL,errFlag] = armaVarCovLS(trajectories,wnVector,length(arParamK(1,:)),length(maParamK(1,:)));
