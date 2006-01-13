function [innovation,innovationVar,wnVector,errFlag] = ...
    armaXKalmanInnov(trajOut,trajIn,arParam,maParam,xParam)
%ARMAXKALMANINNOV finds the innovations (and their variances) resulting from fitting an ARMAX model to a time series using Kalman prediction and filtering.
%
%SYNOPSIS [innovation,innovationVar,wnVector,errFlag] = ...
%    armaXKalmanInnov(trajOut,trajIn,arParam,maParam,xParam)
%
%INPUT  trajOut: Trajectory to be modeled (with measurement uncertainties).
%                Missing points should be indicated with NaN.
%       trajIn : Input series. Missing points should be indicated with NaN.
%                Enter [] if there is no input series.
%       arParam: Autoregressive coefficients (row vector).
%       maParam: Moving average coefficients (row vector).
%       xParam : Coefficients indicating dependence on input (row vector).
%
%OUTPUT innovation   : Vector of differences between predicted and observed data, or innovations.
%       innovationVar: Vector of innovation variances.
%       wnVector     : Estimated white noise in the process.
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

innovaton = [];
innovationVar = [];
wnVector = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if correct number of arguments was used when function was called
if nargin < nargin('armaXKalmanInnov')
    disp('--armaXKalmanInnov: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%find trajectory length
trajLength = size(trajOut,1);

%make sure that input and output series have the same length
if ~isempty(trajIn)
    if size(trajIn,1)~=trajLength
        disp('--armaXKalmanInnov: Input and output series must have the same length!');
        errFlag  = 1;
        return
    end
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

%add maxOrder zeros to end of input series
trajIn = [trajIn; zeros(maxOrder,2)];

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

%construct row vector H (Eq. 2.17)
observationVec = zeros(1,maxOrder);
observationVec(1) = 1;

%initialize state vector and its covariance matrix
stateVecT_T = zeros(maxOrder,1); %Z(0|0)
[stateCovMatT_T,errFlag] = covKalmanInit(arParam,maParam,procErrCov,...
    arOrder,maOrder,maxOrder); %P(0|0)

%initialize innovations vector, its covariance matrix and white noise vector
innovation = NaN*ones(trajLength,1);
innovationVar = NaN*ones(trajLength,1);
wnVector = NaN*ones(trajLength,1);

%initialize vector of contributions from input
inputContr = zeros(maxOrder,1);

%go over all points in trajectory
%note that in the iterations t+1 = timePoint, t = timePoint-1
for timePoint = 1:trajLength

    %calculate contribution of input series (See my notes for derivation)
    if ~isempty(xOrder)
        modXOrder = timePoint - 1 + maxOrder - max(1,timePoint-1+maxOrder-xOrder);
        inputContr(end) = xParam(1:modXOrder+1)*...
            trajIn(timePoint-1+maxOrder:-1:timePoint-1+maxOrder-modXOrder,1);
    end
        
    %predict state at time t+1 given state at time t
    stateVecT1_T = transitionMat*stateVecT_T + inputContr; %Z(t+1|t), Eq. 3.1
    
    %obtain the predicted state's covariance matrix
    stateCovMatT1_T = transitionMat*stateCovMatT_T*transitionMat' ...
        + procErrCov*procErrCov'; %P(t+1|t), Eq. 3.2
    
    %make sure that matrix is symmetric
    stateCovMatT1_T = (stateCovMatT1_T+stateCovMatT1_T')/2;
    
    %predict observable at time t+1
    observableP = stateVecT1_T(1); %y(t+1|t), Eq. 3.3
    
    if isnan(trajOut(timePoint,1)) %if observation at this time point is missing
        
        %cannot modify state vector and its covariance matrix predicted 
        %from previous timepoint since there is no observation
        stateVecT_T = stateVecT1_T; %Z(t+1|t+1), Eq. 5.1
        stateCovMatT_T = stateCovMatT1_T; %P(t+1|t+1), Eq. 5.2
        
    else %if there is an observation
        
        %get innovation, dy(t+1) (Eq. 3.8)
        innovation(timePoint) = trajOut(timePoint,1) - observableP; 
        
        %and its variance, V(t+1) (Eq. 3.6 & 3.10)
        innovationVar(timePoint) = stateCovMatT1_T(1,1) + trajOut(timePoint,2)^2;

        %calculate delta
        delta = stateCovMatT1_T*observationVec'/innovationVar(timePoint); %delta(t+1), Eq. 3.5

        %modify state vector prediction using observation
        stateVecT_T = stateVecT1_T + delta*innovation(timePoint); %Z(t+1|t+1), Eq. 3.4
        
        %update state covariance matrix
        stateCovMatT_T = stateCovMatT1_T - delta*observationVec*stateCovMatT1_T; %P(t+1|t+1), Eq. 3.7
        
        %make sure that matrix is symmetric
        stateCovMatT_T = (stateCovMatT_T+stateCovMatT_T')/2;
        
        %calculate white noise at this time point using 
        %wn(t+1) = x(t+1|t+1) - x(t+1|t) (Eq. 2.10 with j=1)
        wnVector(timePoint) = stateVecT_T(1) - stateVecT1_T(1);
        
    end %(if isnan(trajOut(timePoint,1)) ... else ...)
    
end %(for timePoint = 1:trajLength)


%%%%% ~~ the end ~~ %%%%%
