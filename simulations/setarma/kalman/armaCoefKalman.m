function [arParam,maParam,wnVariance,wnVector,aic,aicc,errFlag] = ...
    armaCoefKalman(trajectories,arParamP0,maParamP0)
%ARMACOEFKALMAN fits an ARMA(p,q) model to a time series which could have missing data points using Kalman prediction and filtering.
%
%SYNOPSIS [arParam,maParam,wnVariance,wnVector,aic,aicc,errFlag] = ...
%    armaCoefKalman(trajectories,arParamP0,maParamP0)
%
%INPUT  trajectories: Structure array containing trajectories to be modeled:
%           .observations: 2D array of measurements and their uncertainties.
%                          Missing points should be indicated with NaN.
%       arParamP0   : Initial guess of partial autoregressive coefficients (row vector).
%       maParamP0   : Initial guess of partial moving average coefficients (row vector).
%
%OUTPUT arParam     : Estimated AR parameters.
%       maParam     : Estimated MA parameters.
%       wnVariance  : Estimated variance of white noise.
%       wnVector    : Structure containing 1 field:
%           .series     : Estimated white noise series in a trajectory.
%       aic         : Akaike's Information Criterion.
%       aicc        : Akaike's Information Criterion with bias Correction.
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%REMARKS The algorithm implemented here is that presented in R. H. Jones,
%        "Maximum Likelihood Fitting of ARMA Models to Time Series with
%        Missing Observations", Technometrics 22: 389-395 (1980). All
%        equation numbers used here are those in that paper. The main
%        difference is that I do not estimate the observational error
%        variance, but use that obtained from experimental data or
%        simulated trajectories (and is thus time-dependent). Also, 
%        trajectories are shifted to get zero mean before analysis is done.
%
%Khuloud Jaqaman, July 2004

%initialize output
arParam = [];
maParam = [];
wnVariance = [];
wnVector = [];
aic = [];
aicc = [];
errFlag = 0;

%check if correct number of arguments was used when function was called
if nargin < nargin('armaCoefKalman')
    disp('--armaCoefKalman: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check input data
for i=1:length(trajectories);
    traj = trajectories(i).observations;
    [trajLength,nCol] = size(traj);
    if nCol ~= 2
        if nCol == 1 %if no error is supplied, it is assumed that there is no observational error
            traj = [traj zeros(trajLength,1)];
        else
            disp('--armaCoefKalman: "trajectories.observations" should have either 1 column for measurements, or 2 columns: 1 for measurements and 1 for measurement uncertainties!');
            errFlag = 1;
        end
    end
    trajectories(i).observations = traj;
end
[nRow,arOrder] = size(arParamP0);
if ~isempty(arParamP0)
    if nRow ~= 1
        disp('--armaCoefKalman: "arParamP0" should be a row vector!');
        errFlag = 1;
    end
    if ~isempty(find(abs(arParamP0)>=1))
        disp('--armaCoefKalman: All entries in "arParamP0" should be smaller than 1 in magnitude!');
        errFlag = 1;
    end
end
[nRow,maOrder] = size(maParamP0);
if ~isempty(maParamP0)
    if nRow ~= 1
        disp('--armaCoefKalman: "maParamP0" should be a row vector!');
        errFlag = 1;
    end
    if ~isempty(find(abs(maParamP0)>=1))
        disp('--armaCoefKalman: All entries in "maParamP0" should be smaller than 1 in magnitude!');
        errFlag = 1;
    end
end
if errFlag
    disp('--armaCoefKalman: Please fix input data!');
    return
end

%initial values of parameters to be estimated
param0 = [arParamP0 maParamP0];

%obtain number of available observations and their indices, and shift
%trajectory to get zero mean
for i=1:length(trajectories)
    traj = trajectories(i).observations;
    trajectories(i).available = find(~isnan(traj(:,1))); %get indices of available points
    numAvail(i) = length(trajectories(i).available); %get total # of available points
    traj(:,1) = traj(:,1) - mean(traj(trajectories(i).available,1)); %shift trajectory
    trajectories(i).observations = traj;
end
%calculcate overall number of available points (in all trajectories)
totAvail = sum(numAvail); 

%define optimization options.
options = optimset('Display','final','DiffMaxChange',1e-3,...
    'DiffMinChange',1e-8,'TolFun',1e-4,'TolX',1e-4); 

%minimize -2ln(likelihood) to get ARMA coefficients and variance of observational error
[params,fval,exitFlag] = fmincon(@armaCoefKalmanObj,param0,[],[],[],[],...
    -0.99*ones(1,arOrder+maOrder),0.99*ones(1,arOrder+maOrder),...
    [],options,arOrder,trajectories);

%if minimization was successful
if exitFlag > 0

    %assign parameters obtained through minimization
    arParamP = params(1:arOrder);
    maParamP = params(arOrder+1:end);

    %get AR and MA coefficients
    if ~isempty(arParamP)
        [arParam,errFlag] = levinsonDurbinAR(arParamP);
    else
        arParam = [];
    end
    if ~isempty(maParamP)
        [maParam,errFlag] = levinsonDurbinMA(maParamP);
    else
        maParam = [];
    end

    %check for causality and invertibility of estimated model
    r = abs(roots([-arParam(end:-1:1) 1]));
    if ~isempty(find(r<=1.00001))
        disp('--armaCoefKalman: Predicted model not causal!');
        errFlag = 1;
        return
    end
    r = abs(roots([maParam(end:-1:1) 1]));
    if ~isempty(find(r<=1.00001))
        disp('--armaCoefKalman: Predicted model not invertible!');
        errFlag = 1;
        return
    end

    %obtain likelihood, white noise sequence and white noise variance
    sum1 = 0;
    sum2 = 0;
    for i = 1:length(trajectories)

        %obtain available points in this trajectory
        available = trajectories(i).available;

        %get the innovations, their variances and the estimated white noise series
        %using Kalman prediction and filtering
        [innovation,innovationVar,wnVector(i).series,errFlag] = ...
            armaKalmanInnov(trajectories(i).observations,arParam,maParam);
        if errFlag
            disp('--armaCoefKalman: "armaKalmanInnov" did not function properly!');
            return
        end

        %calculate white noise variance of current trajectory
        wnVarianceSamp(i) = mean(innovation(available).^2./innovationVar(available));

        %1st sum in Eq. 3.15
        sum1 = sum1 + sum(log(innovationVar(available)));
        %2nd sum in Eq. 3.15
        sum2 = sum2 + sum(innovation(available).^2./innovationVar(available));

    end %(for i = 1:length(trajectories))

    %calculate -2ln(likelihood)
    neg2LnLikelihood = sum1 + totAvail*log(sum2);

    %calculate mean white noise variance of all trajectories
    wnVariance = sum(numAvail.*wnVarianceSamp)/totAvail;

    %get number of parameters estimated: arOrder AR coefficients, maOrder MA
    %coefficients and white noise variance
    numParam = arOrder + maOrder + 1;

    %evaluate Akaike's Information Criterion
    aic = neg2LnLikelihood + 2*numParam;

    %evaluate the bias-corrected Akaike's Information Criterion
    aicc = neg2LnLikelihood + 2*numParam*totAvail/(totAvail-numParam-1);

    %put partial coefficients in 2nd row of coefficients matrix
    arParam(2,:) = arParamP;
    maParam(2,:) = maParamP;

else %if minimization was not successful
    
    errFlag = 1;
    
end %(if exitFlag > 0)
    
% % %TOMLAB stuff
% % % prob = glcAssign('armaCoefKalmanObj',-0.99*ones(1,arOrder+maOrder),...
% % %     0.99*ones(1,arOrder+maOrder),'minNegLik',[],[],[],[],[],[],param0);
% % % prob.user.arOrder = arOrder;
% % % prob.user.trajectories = trajectories;
% % % prob.PriLevOpt = 15;
% % % prob.optParam.MaxFunc = 3000;
% % % prob.optParam.MaxIter = 3000;
% % % result = tomRun('glbFast',prob,0,2);
% % % result2 = tomRun('conSolve',prob,0,2);
% % % params = result.x_k(:,1)';
% % % params2 = result.x_k(:,1)';


