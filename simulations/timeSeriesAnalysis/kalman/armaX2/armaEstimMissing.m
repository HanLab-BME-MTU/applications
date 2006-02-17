function [trajFilled,errFlag] = armaEstimMissing(trajOut,arParam,maParam,...
    wnVariance)
%ARMAESTIMMISSING estimates missing observations in an ARMA series using fixed-interval smoothing
%
%SYNOPSIS [trajFilled,errFlag] = armaEstimMissing(trajOut,arParam,maParam,...
%    wnVariance)
%
%INPUT  trajOut   : Observations of output time series to be fitted. Either an
%                   array of structures trajOut(1:nTraj).observations, or a
%                   2D array representing one single trajectory.
%           .observations: 2D array of measurements and their uncertainties.
%                   Missing points should be indicated with NaN.
%       arParam   : Autoregressive coefficients (row vector). Enter as []
%                   if there is no AR part.
%       maParam   : Moving average coefficients (row vector). Enter as []
%                   if there is no MA part.
%       wnVariance: White noise variance.
%
%OUTPUT trajFilled: Trajectories with estimated values of missing
%                   observations. Structure with fields observations.
%       errFlag   : 0 if function executes normally, 1 otherwise.
%
%REMARKS Smoothing algorithm taken from Andrew C. Harvey, "Forecasting,
%        structural time series model and the Kalman filter," Section 3.6.2
%        (p. 154). State-space formulation of an ARMA model taken from
%        Jones's paper.
%
%MATLAB VERSION (originally written on): 7.0.4.352 (R14) Service Pack 2 
%
%Khuloud Jaqaman, February 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trajFilled = [];
errFlag = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if correct number of arguments was used when function was called
if nargin < nargin('armaEstimMissing')
    disp('--armaEstimMissing: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check trajOut and turn it into struct if necessary
if ~isstruct(trajOut)
    tmp = trajOut;
    clear trajOut
    trajOut.observations = tmp;
    clear tmp
elseif ~isfield(trajOut,'observations')
    disp('--armaEstimMissing: Please input "trajOut" in fields ''observations''!')
    errFlag = 1;
    return
end

%get number of trajectories
numTraj = length(trajOut);

%make a copy of trajOut in trajFilled
trajFilled = trajOut;

%get length of each trajectory and add column of zeros for observational 
%error if not input
trajLength = zeros(numTraj,1);
for i=1:numTraj

    trajT = trajOut(i).observations;
    [trajLength(i),nCol] = size(trajT);
    
    if nCol ~= 2
        if nCol == 1 %if no error is supplied, it is assumed that there is no observational error
            trajT = [trajT zeros(trajLength,1)];
        else
            disp('--armaEstimMissing: "traj.observations" should have either 1 column for measurements, or 2 columns: 1 for measurements and 1 for measurement uncertainties!');
            errFlag = 1;
        end
    end
    trajOut(i).observations = trajT;

end

%get arOrder
[nRow,arOrder] = size(arParam);
if ~isempty(arParam)
    if nRow ~= 1
        disp('--armaEstimMissing: "arParam" should be a row vector!');
        errFlag = 1;
    end
end

%get maOrder
[nRow,maOrder] = size(maParam);
if ~isempty(maParam)
    if nRow ~= 1
        disp('--armaEstimMissing: "maParam" should be a row vector!');
        errFlag = 1;
    end
end

%check white noise variance
if wnVariance < 0
    disp('--armaEstimMissing: "wnVariance" should be nonnegative!');
    errFlag = 1;
end

%exit if there are problems in input data
if errFlag
    disp('--armaEstimMissing: Please fix input data!');
    return
end

%get maxOrder to determine size of matrices and vectors in Eq. 2.15 - 2.17
maxOrder = max(arOrder,maOrder+1);

%add zeros to ends of arParam and maParam to get vectors of length maxOrder
arParamMod = [arParam zeros(1,maxOrder-arOrder)];
maParamMod = [maParam zeros(1,maxOrder-maOrder)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Smoothing algorithm
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

%get initial state vector and its covariance matrix
stateVec00 = zeros(maxOrder,1); %Z(0|0)
[stateCovMat00,errFlag] = covKalmanInit(arParam,maParam,procErrCov,...
    arOrder,maOrder,maxOrder); %P(0|0)
stateCovMat00 = stateCovMat00*wnVariance;

%go over all supplied trajectories
for itraj=1:numTraj

    %get trajectories
    trajOutI = trajOut(itraj).observations;

    %reserve memory for filtered state and its covariance matrix at all time points
    stateVecP = NaN*ones(maxOrder,trajLength(itraj));
    stateVecF = NaN*ones(maxOrder,trajLength(itraj));
    stateCovP = NaN*ones(maxOrder,maxOrder,trajLength(itraj));
    stateCovF = NaN*ones(maxOrder,maxOrder,trajLength(itraj));

    %initialize state vector and its covariance matrix
    stateVecT_T = stateVec00;
    stateCovMatT_T = stateCovMat00;

    %first get all predicted and filtered estimates for all points in trajectory
    %note that in the iterations t+1 = timePoint, t = timePoint-1
    for timePoint = 1:trajLength(itraj)

        %predict state at time t+1 given state at time t
        stateVecT1_T = transitionMat*stateVecT_T; %Z(t+1|t), Eq. 3.1

        %obtain the predicted state's covariance matrix
        stateCovMatT1_T = transitionMat*stateCovMatT_T*transitionMat' ...
            + wnVarianceGG; %P(t+1|t), Eq. 3.2

        %save predicted state and covariance matrix
        stateVecP(:,timePoint) = stateVecT1_T;
        stateCovP(:,:,timePoint) = stateCovMatT1_T;

        if isnan(trajOutI(timePoint,1)) %if observation at this time point is missing

            %cannot modify state vector and its covariance matrix predicted
            %from previous timepoint since there is no observation
            stateVecT_T = stateVecT1_T; %Z(t+1|t+1), Eq. 5.1
            stateCovMatT_T = stateCovMatT1_T; %P(t+1|t+1), Eq. 5.2

            %make sure that covariance matrix is symmetric
            stateCovMatT_T = (stateCovMatT_T+stateCovMatT_T')/2;

        else %if there is an observation

            %get innovation, dy(t+1) (Eq. 3.8)
            innovation = trajOutI(timePoint,1) - stateVecT1_T(1);

            %and its variance, V(t+1) (Eq. 3.6 & 3.10)
            innovationVar = stateCovMatT1_T(1,1) + trajOutI(timePoint,2)^2;

            %calculate delta
            delta = stateCovMatT1_T(:,1)/innovationVar; %delta(t+1), Eq. 3.5

            %modify state vector prediction using observation
            stateVecT_T = stateVecT1_T + delta*innovation; %Z(t+1|t+1), Eq. 3.4

            %update state covariance matrix
            stateCovMatT_T = stateCovMatT1_T - delta*observationVec*stateCovMatT1_T; %P(t+1|t+1), Eq. 3.7

            %make sure that matrix is symmetric
            stateCovMatT_T = (stateCovMatT_T+stateCovMatT_T')/2;

        end %(if isnan(trajOutI(timePoint,1)) ... else ...)

        %store value of forward-predicted state vector and its covariance matrix
        stateVecF(:,timePoint) = stateVecT_T;
        stateCovF(:,:,timePoint) = stateCovMatT_T;
        
    end %(for timePoint = 1:trajLength(itraj))

    %get filtered state at last time point
    stateVecT_N = stateVecF(:,end);
    stateCovT_N = stateCovF(:,:,end);
    
    %then calculate the best estimate of the missing observations using all 
    %available observations
    for timePoint=trajLength(itraj)-1:-1:1
        
        %calculate intermediate matrix
        interMat = stateCovF(:,:,timePoint)*transitionMat'*...
            inv(stateCovP(:,:,timePoint+1));

        %get the best estimate of the state vector
        stateVecT_N = stateVecF(:,timePoint) + interMat*(stateVecT_N - ...
            transitionMat*stateVecF(:,timePoint));

        %get the covariance matrix
        stateCovT_N = stateCovF(:,:,timePoint) + interMat*(stateCovT_N - ...
            stateCovP(:,:,timePoint+1))*interMat';
            
        %get the best estimate of the missing observation
        if isnan(trajFilled(itraj).observations(timePoint,1))
            trajFilled(itraj).observations(timePoint,1) = ...
                stateVecT_N(1);
        end

    end %(for timePoint=trajLength(itraj)-1:-1:1)

end %(for itraj=1:numTraj)


%%%%% ~~ the end ~~ %%%%%
