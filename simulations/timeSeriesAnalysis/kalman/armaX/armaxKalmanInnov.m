function [innovation,innovationVar,wnVector,stateVec,stateCov,errFlag] = ...
    armaxKalmanInnov(trajOut,trajIn,arParam,maParam,xParam,wnVariance)
%ARMAXKALMANINNOV does forward Kalman prediction and filtering of a time series using an ARMAX model
%
%SYNOPSIS [innovation,innovationVar,wnVector,stateVec,stateCov,errFlag] = ...
%    armaxKalmanInnov(trajOut,trajIn,arParam,maParam,xParam,wnVariance)
%
%INPUT  trajOut   : Trajectory to be modeled (with measurement uncertainties).
%                   Missing points should be indicated with NaN.
%       trajIn    : Input series. Must not have any missing points.
%                   Enter [] if there is no input series.
%       arParam   : Autoregressive coefficients (row vector).
%       maParam   : Moving average coefficients (row vector).
%       xParam    : Coefficients indicating dependence on input (row vector).
%       wnVariance: White noise Variance. Optional. Default: 1.
%
%OUTPUT innovation   : Vector of differences between predicted and observed data, or innovations.
%       innovationVar: Vector of innovation variances.
%       wnVector     : Estimated white noise in the process.
%       stateVec     : Forward-predicted state vector at missing time points.
%       stateCov     : Forward-predicted covariance matrix at missing time points.
%       errFlag      : 0 if function executes normally, 1 otherwise.
%
%REMARKS The algorithm implemented here is an ARMAX generalized version of
%        the algorithm presented in R. H. Jones,"Maximum Likelihood Fitting 
%        of ARMA Models to Time Series with Missing Observations", 
%        Technometrics 22: 389-395 (1980). All equation numbers used are 
%        those in the paper. However, I do not estimate the observational 
%        error variance, but use that obtained from experimental data or
%        simulated trajectories (and is thus time-dependent).
%
%Khuloud Jaqaman, January 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

innovation = [];
innovationVar = [];
wnVector = [];
stateVec = [];
stateCov = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if correct number of arguments was used when function was called
if nargin < 5
    disp('--armaxKalmanInnov: Incorrect number of input arguments!');
    errFlag = 1;
    return
end

%find trajectory length and number of missing observations
trajLength = size(trajOut,1);
numMissing = length(find(isnan(trajOut(:,1))));

%make sure that input and output series have the same length
if ~isempty(trajIn)
    if size(trajIn,1)~=trajLength
        disp('--armaxKalmanInnov: Input and output series must have the same length!');
        errFlag = 1;
        return
    end
end

%assign 1 to WN variance if not supplied
if nargin < 6 || isempty(wnVariance)
    wnVariance = 1;
end

%find arOrder, maOrder and xOrder
arOrder = length(arParam);
maOrder = length(maParam);
xOrder  = length(xParam) - 1;

%get maxOrder to determine size of matrices and vectors in Eq. 2.15 - 2.17
maxOrder = max(arOrder,maOrder+1);

%add zeros to ends of arParam and maParam to get vectors of length maxOrder
arParamMod = [arParam zeros(1,maxOrder-arOrder)];
maParamMod = [maParam zeros(1,maxOrder-maOrder)];

%rewrite dependence on input
if xOrder == -1 %if there is no dependence
    xOrder = 0;
    xParam = 0;
    trajIn = zeros(trajLength+maxOrder,2);
else %if there is dependence
    trajIn = [zeros(xOrder,2); trajIn; zeros(maxOrder,2)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculation of innovations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%construct matrix F (Eqs. 2.15, 2.16)
transitionMat = diag(ones(maxOrder-1,1),1);
transitionMat(end,:) = arParamMod(end:-1:1); 

%construct column vector G (Eqs. 2.15, 2.16) using the recursions in Eq. 2.13
%note that procErrCov(i) is the covariance of the process at time t 
%with the white noise at time t+i-1, normalized by the white noise variance
%(Eqs. 4.2, 4.3 and 4.4)
procErrCov = ones(maxOrder,1);
for i = 2:maxOrder
    dummy = maParamMod(i-1) + arParamMod(1:i-1)*procErrCov(i-1:-1:1);
    procErrCov(i) = dummy;
end

%calculate wnVariance*G*G'
wnVarianceGG = wnVariance*procErrCov*procErrCov';

%construct row vector H (Eq. 2.17)
observationVec = zeros(1,maxOrder);
observationVec(1) = 1;

%construct matrix of dependence on input
inputCoefMat = zeros(maxOrder,xOrder+1);
inputCoefMat(end,:) = xParam(end:-1:1);

%initialize state vector and its covariance matrix
stateVecT_T = zeros(maxOrder,1); %Z(0|0)
[stateCovMatT_T,errFlag] = covKalmanInit(arParam,maParam,procErrCov,...
    arOrder,maOrder,maxOrder); %P(0|0)
stateCovMatT_T = stateCovMatT_T*wnVariance;

%initialize output variables and indxMiss
innovation = NaN*ones(trajLength,1);
innovationVar = NaN*ones(trajLength,1);
wnVector = NaN*ones(trajLength,1);
% stateVec = NaN*ones(maxOrder,numMissing);
% stateCov = NaN*ones(maxOrder,maxOrder,numMissing);
stateVec = NaN*ones(maxOrder,trajLength);
stateCov = NaN*ones(maxOrder,maxOrder,trajLength);
indxMiss = 0;

%go over all points in trajectory
%note that in the iterations t+1 = timePoint, t = timePoint-1
% % % % for timePoint = 1:trajLength
for timePoint = max(xOrder+1,1):trajLength

    %get vector of inputs
    tmp = timePoint - 1 + maxOrder;
    inputVector = trajIn(tmp:tmp+xOrder,1);

    %predict state at time t+1 given state at time t
    stateVecT1_T = transitionMat*stateVecT_T + inputCoefMat*inputVector; %Z(t+1|t), Eq. 3.1
   
    %obtain the predicted state's covariance matrix
    stateCovMatT1_T = transitionMat*stateCovMatT_T*transitionMat' ...
        + wnVarianceGG; %P(t+1|t), Eq. 3.2
    
    if isnan(trajOut(timePoint,1)) %if observation at this time point is missing
        
        %cannot modify state vector and its covariance matrix predicted 
        %from previous timepoint since there is no observation
        stateVecT_T = stateVecT1_T; %Z(t+1|t+1), Eq. 5.1
        stateCovMatT_T = stateCovMatT1_T; %P(t+1|t+1), Eq. 5.2

        %make sure that covariance matrix is symmetric
        stateCovMatT_T = (stateCovMatT_T+stateCovMatT_T')/2;
        
        %store values of forward-predicted state vector and its covariance matrix
%         indxMiss = indxMiss + 1;
%         stateVec(:,indxMiss) = stateVecT_T;
%         stateCov(:,:,indxMiss) = stateCovMatT_T;
        stateVec(:,timePoint) = stateVecT_T;
        stateCov(:,:,timePoint) = stateCovMatT_T;
        
    else %if there is an observation

        %get innovation, dy(t+1) (Eq. 3.8)
        innovation(timePoint) = trajOut(timePoint,1) - stateVecT1_T(1);

        %and its variance, V(t+1) (Eq. 3.6 & 3.10)
        innovationVar(timePoint) = stateCovMatT1_T(1,1) + trajOut(timePoint,2)^2;

        %calculate delta
        delta = stateCovMatT1_T(:,1)/innovationVar(timePoint); %delta(t+1), Eq. 3.5

        %modify state vector prediction using observation
        stateVecT_T = stateVecT1_T + delta*innovation(timePoint); %Z(t+1|t+1), Eq. 3.4
        
        %update state covariance matrix
        stateCovMatT_T = stateCovMatT1_T - delta*observationVec*stateCovMatT1_T; %P(t+1|t+1), Eq. 3.7
        
        %make sure that matrix is symmetric
        stateCovMatT_T = (stateCovMatT_T+stateCovMatT_T')/2;
        
        %calculate white noise at this time point using 
        %wn(t+1) = x(t+1|t+1) - x(t+1|t) (Eq. 2.10 with j=1)
        wnVector(timePoint) = stateVecT_T(1) - stateVecT1_T(1);
        
        %store values of forward-predicted state vector and its covariance matrix
        stateVec(:,timePoint) = stateVecT_T;
        stateCov(:,:,timePoint) = stateCovMatT_T;

    end %(if isnan(trajOut(timePoint,1)) ... else ...)
    
end %(for timePoint = 1:trajLength)


%%%%% ~~ the end ~~ %%%%%
