function [selectCrit,errFlag] = selectCritKalman(trajectories,arParam,...
    maParam)
%SELECTCRITKALMAN calculates the AIC, AICC and BIC of a model wrt a set of trajectories
%
%SYNOPSIS [selectCrit,errFlag] = selectCritKalman(trajectories,arParam,...
%    maParam)
%
%INPUT  trajectories: Observations of time series to be fitted. Either an 
%                     array of structures traj(1:nTraj).observations, or a
%                     2D array representing one single trajectory. 
%           .observations: 2D array of measurements and their uncertainties.
%                     Missing points should be indicated with NaN.
%       arParam     : Autoregressive coefficients in model (row vector).
%       maParam     : Moving average coefficients in model (row vector).
%
%OUTPUT selectCrit  : Model selection criterion. Contains 3 fields:
%            .aic        : Akaike's Information Criterion.
%            .aicc       : Bias-Corrected Akaike's Information Criterion.
%            .bic        : Bayesian Information Criterion.
%
%REMARKS All equation numbers used here refer to R. H. Jones, "Maximum 
%        Likelihood Fitting of ARMA Models to Time Series with Missing
%        Observations", Technometrics 22: 389-395 (1980). 
%
%Khuloud Jaqaman, October 2004

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

selectCrit = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin ~= 3
    disp('--selectCritKalman: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check trajectory and turn it into struct if necessary
if ~isstruct(trajectories)
    tmp = trajectories;
    clear trajectories
    trajectories.observations = tmp;
    clear tmp
elseif ~isfield(trajectories,'observations')
    disp('--autoCorr: Please input the trajectories in fields ''observations''')
    errFlag = 1;
    return
end

%get number of trajectories and add column for observational error if necessary
trajOriginal = trajectories;
for i=1:length(trajectories);
    traj = trajectories(i).observations;
    [trajLength,nCol] = size(traj);
    if nCol ~= 2
        if nCol == 1 %if no error is supplied, it is assumed that there is no observational error
            traj = [traj zeros(trajLength,1)];
        else
            disp('--selectCritKalman: "trajectories.observations" should have either 1 column for measurements, or 2 columns: 1 for measurements and 1 for measurement uncertainties!');
            errFlag = 1;
        end
    end
    trajectories(i).observations = traj;
end

%get number of AR coefficients
[nRow,arOrder] = size(arParam);
if ~isempty(arParam)
    if nRow ~= 1
        disp('--selectCritKalman: "arParam" should be a row vector!');
        errFlag = 1;
    end
end

%get number of MA coefficients
[nRow,maOrder] = size(maParam);
if ~isempty(maParam)
    if nRow ~= 1
        disp('--selectCritKalman: "maParam" should be a row vector!');
        errFlag = 1;
    end
end

%get number of parameters "estimated"
numParam = arOrder + maOrder + 1;

%exit function if there are problems in input data
if errFlag
    disp('--selectCritKalman: Please fix input data!');
    return
end

%obtain number of available observations and shift trajectories to get zero mean
for i=1:length(trajectories)
    
    traj = trajectories(i).observations;
    numAvail(i) = length(find(~isnan(traj(:,1)))); %get # of available points
    traj(:,1) = traj(:,1) - nanmean(traj(:,1)); %shift trajectory
    trajectories(i).observations = traj;

end
totAvail = sum(numAvail); %calculate total number of available points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Likelihood calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%go over all trajectories and calculate overall likelihood
sum1 = 0;
sum2 = 0;
for i = 1:length(trajectories)

    %get the innovations, their variances and the estimated white noise series
    %using Kalman prediction and filtering
    [innovation,innovationVar,wnVector(i).observations,dummy1,dummy2,errFlag] = ...
        armaKalmanInnov(trajectories(i).observations,arParam,maParam);
    if errFlag
        disp('--selectCritKalman: "armaKalmanInnov" did not function properly!');
        return
    end

    %1st sum in Eq. 3.15
    sum1 = sum1 + nansum(log(innovationVar));
    %2nd sum in Eq. 3.15
    sum2 = sum2 + nansum(innovation.^2./innovationVar);

end %(for i = 1:length(trajectories))

%calculate -2ln(likelihood)
neg2LnLikelihoodV = sum1 + totAvail*log(sum2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Selection criteria evaluation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Akaike's Information Criterion
selectCrit.aic = neg2LnLikelihoodV + 2*numParam;

%Bias-corrected Akaike's Information Criterion
selectCrit.aicc = neg2LnLikelihoodV + 2*numParam*totAvail/(totAvail-numParam-1);

%Bayesian Information Criterion
selectCrit.bic = neg2LnLikelihoodV + log(totAvail)*numParam;
    

%%%%% ~~ the end ~~ %%%%%
