function [innovation,innovationVar,wnVector,errFlag] = ...
    armaKalmanInnov(traj,arParam,maParam);
%ARMAKALMANINNOV finds the innovations (and their variances) resulting from fitting an ARMA(p,q) model to a time series which could have missing data points using Kalman prediction and filtering.
%
%SYNOPSIS [innovation,innovationVar,wnVector,errFlag] = ...
%    armaKalmanInnov(traj,arParam,maParam);
%
%INPUT  traj       : Trajectory to be modeled (with measurement uncertainties).
%                    Missing points should be indicated with NaN.
%       arParam    : Autoregressive coefficients (row vector).
%       maParam    : Moving average coefficients (row vector).
%
%OUTPUT innovation   : Vector of differences between predicted and observed data, or innovations.
%       innovationVar: Vector of innovation variances.
%       wnVector     : Estimated white noise in the process.
%       errFlag      : 0 if function executes normally, 1 otherwise.
%
%REMARKS The algorithm implemented here is that presented in R. H. Jones,
%        "Maximum Likelihood Fitting of ARMA Models to Time Series with
%        Missing Observations", Technometrics 22: 389-395 (1980). All
%        equation numbers used here are those in that paper. The main
%        difference is that I do not estimate the observational error
%        variance, but use that obtained from experimental data or
%        simulated trajectories (and is thus time-dependent).
%
%Khuloud Jaqaman, July 2004

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

innovaton = [];
innovationVar = [];
wnVector = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if correct number of arguments was used when function was called
if nargin < nargin('armaKalmanInnov')
    disp('--armaKalmanInnov: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%find trajectory length, arOrder and maOrder
trajLength = size(traj,1);
arOrder = length(arParam);
maOrder = length(maParam);

%get maxOrder to determine size of matrices in Eq. 2.15 - 2.17
maxOrder = max(arOrder,maOrder+1);

%add zeros to ends of arParam and maParam to get vectors of length maxOrder
arParamMod = [arParam zeros(1,maxOrder-arOrder)];
maParamMod = [maParam zeros(1,maxOrder-maOrder)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculation of innovations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%construct matrix F (Eq. 2.15, 2.16)
transitionMat = diag(ones(maxOrder-1,1),1);
transitionMat(end,:) = arParamMod(end:-1:1); 

%construct column vector G (Eq. 2.15, 2.16) using the recursions in Eq. 2.13
%note that procErrCov(i) is the covariance of the process at time t 
%with the white noise at time t+i-1, normalized by the white noise variance
%(Eqs. 4.2, 4.3 and 4.4)
procErrCov = ones(maxOrder,1);
for i = 2:maxOrder
    dummy = maParamMod(i-1) + arParamMod(1:i-1)*procErrCov(i-1:-1:1);
    procErrCov(i) = dummy;
end

%construct row vector H (Eq. 2.17)
observationVec = zeros(1,maxOrder);
observationVec(1) = 1;

%initialize state vector and its covariance matrix
stateVecT_T = zeros(maxOrder,1); %Z(0|0)
% [stateCovMatT_T,errFlag] = covKalmanInit(arParam,maParam,procErrCov,...
%     arOrder,maOrder,maxOrder); %P(0|0)
stateCovMatT_T = 1e4*ones(maxOrder);

%initialize innovations vector, its covariance matrix and white noise vector
innovation = NaN*ones(trajLength,1);
innovationVar = NaN*ones(trajLength,1);
wnVector = NaN*ones(trajLength,1);

%got over all points in trajectory
for timePoint = 1:trajLength
    
    %predict state at time T+1 given state at time T
    stateVecT1_T = transitionMat*stateVecT_T; %Z(t+1|t), Eq. 3.1
    %obtain the predicted state's covariance matrix
    stateCovMatT1_T = transitionMat*stateCovMatT_T*transitionMat' ...
        + procErrCov*procErrCov'; %P(t+1|t), Eq. 3.2
    %predict observable at time T+1
    observableP = stateVecT1_T(1); %y(t+1|t), Eq. 3.3
    
    if isnan(traj(timePoint,1)) %if observation at this time point is missing
        
        %cannot modify state vector and its covariance matrix predicted 
        %from previous timepoint since there is no observation
        stateVecT_T = stateVecT1_T; %Z(t+1|t+1), Eq. 5.1
        stateCovMatT_T = stateCovMatT1_T; %P(t+1|t+1), Eq. 5.2
        
    else %if there is an observation
        
        %get innovation, dy(t+1) (Eq. 3.8)
        innovation(timePoint) = traj(timePoint,1) - observableP; 
        %and its variance, V(t+1) (Eq. 3.6 & 3.10)
        innovationVar(timePoint) = stateCovMatT1_T(1,1) + traj(timePoint,2)^2;

        if innovationVar(timePoint) < 0
            disp('neg!');
%             errFlag = 1;
%             return
        end
        
        %calculate delta
        delta = stateCovMatT1_T*observationVec'/innovationVar(timePoint); %delta(t+1), Eq. 3.5
        %modify state vector prediction using observation
        stateVecT_T = stateVecT1_T + delta*innovation(timePoint); %Z(t+1|t+1), Eq. 3.4
        %update state covariance matrix
        stateCovMatT_T = stateCovMatT1_T - delta*observationVec*stateCovMatT1_T; %P(t+1|t+1), Eq. 3.7
        
        %calculate white noise at this time point using 
        %wn(t+1) = x(t+1|t+1) - x(t+1|t) (Eq. 2.10 with j=1)
        wnVector(timePoint) = stateVecT_T(1) - stateVecT1_T(1);
        
    end %(if isnan(traj(timePoint,1)))
    
end %(for timePoint = 1:trajLength)


%%%%% ~~ the end ~~ %%%%%
