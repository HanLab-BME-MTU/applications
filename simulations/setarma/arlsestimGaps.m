function [arParam0,varCovMat0,noiseSigma0,arParam,varCovMat,noiseSigma,trajP,errFlag] = ...
    arlsestimGaps(traj,arOrder,arTol)
%ARLSESTIMGAPS fits an AR model to a trajectory with missing data points using least square fitting
%
%SYNOPSIS [arParam0,varCovMat0,noiseSigma0,arParam,varCovMat,noiseSigma,trajP,errFlag] = ...
%   arlsestimGaps(traj,arOrder,arTol)
%
%INPUT  traj    : Trajectory to be modeled (with measurement uncertainty).
%                 Missing points should be indicated with NaN.
%       arOrder : Order of proposed AR model.
%       arTol   : Change in AR coefficients below which it is safe to stop.
%
%OUTPUT arParam0    : Parameters estimated by omitting missing data points
%                     and points depending on them.
%       varCovMat0  : Variance-covariance matrix of arParam0.
%       noiseSigma0 : Estimated standard deviation of white noise using AR model with
%                     parameters arParam0 and only points used in their estimation.
%       arParam     : Parameters estimated by filling in the gaps and then using
%                     all points except for the missing ones.
%       varCovMat   : Variance-covariance matrix of arParam.
%       noiseSigma  : Estimated standard deviation of white noise using AR
%                     model with parameters arParam and "non-missing" points.
%       trajP       : Trajectory with estimates of missing points and their uncertainties.
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, March 2004

errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= nargin('arlsestimGaps')
    disp('--arlsestimGaps: Incorrect number of input arguments!');
    errFlag  = 1;
    arParam0 = [];
    varCovMat0 = [];
    noiseSigma0 = [];
    arParam = [];
    varCovMat = [];
    noiseSigma = [];
    trajP = []
    return
end

%check input data
if arOrder < 1
    disp('--arlsestimGaps: Variable "arOrder" should be >= 1!');
    errFlag = 1;
end
[trajLength,nCol] = size(traj);
if trajLength < arOrder
    disp('--arlsestimGaps: Length of trajectory should be larger than model order!');
    errFlag = 1;
end
if nCol ~= 2
    if nCol == 1
        traj = [traj ones(trajLength,1)];
    else
        disp('--arlsestimGaps: "traj" should have either 1 column for measurements, or 2 columns: 1 for measurements and 1 for measurement uncertainties!');
        errFlag = 1;
    end
end
if arTol <= 0
    disp('--arlsestimGaps: Variable "arTol" should be positive!');
    errFlag = 1;
end

%exit if there are problems in input data
if errFlag
    disp('--arlsestimGaps: please fix input data!');
    arParam0 = [];
    varCovMat0 = [];
    noiseSigma0 = [];
    arParam = [];
    varCovMat = [];
    noiseSigma = [];
    trajP = [];
    return
end

%assign zeros to missing data points and ones to existing measurements.
available = ~isnan(traj(:,1));

%get indices of missing points
indx = find(~available);

%get indices of available points
indxP = find(available); 

%check indx - algorithm cannot proceed if there are missing points for time < arOrder.
if ~isempty(indx)
    if indx(1) <= arOrder
        disp('--arlsestimGaps: There are missing points for time points < arOrder!');
        disp('  Please fix input data such that at least the first arOrder time points are available!');
        errFlag = 1;
        arParam0 = [];
        varCovMat0 = [];
        noiseSigma0 = [];
        arParam = [];
        varCovMat = [];
        noiseSigma = [];
        trajP = [];
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate arParam0 using available data points that have enough previous points to
%predict them using the proposed order of the AR model. Note that if there
%are no missing data points, this is basically the solution!

%find data points to be used in fitting
fitSet = indxP(find((indxP(arOrder+1:end)-indxP(1:end-arOrder))==arOrder)+arOrder);
fitLength = length(fitSet);

%construct weighted matrix of previous points multiplying AR coefficients (on LHS of equation)
%[size: fitLength by arOrder]
lhsMat = zeros(fitLength,arOrder);
for i=1:arOrder
    lhsMat(:,i) = traj(fitSet-i,1)./traj(fitSet,2);
end

%construct weighted vector of points used to determine parameters (on RHS of equation)
%[size: fitLength by 1]
rhsVec = traj(fitSet,1)./traj(fitSet,2);

%calculate arParam0 by least square fitting
arParam0 = (lhsMat\rhsVec)';

%check for causality of estimated model
r = abs(roots([-arParam0(end:-1:1) 1]));
if ~isempty(find(r<=1.00001))
    disp('--arlsestimGaps: Warning: Predicted model (arParam0) not causal!');
end

%get vector of weighted residuals
epsilon = lhsMat*arParam0' - rhsVec;

%calculate variance-covariance matrix
varCovMat0 = inv(lhsMat'*lhsMat)*(epsilon'*epsilon/(fitLength-arOrder));

%get standard deviation of white noise (note that nonweighted residuals
%must be used in this calculation)
noiseSigma0 = std(epsilon.*traj(fitSet,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if there are missing data points, predict them using arParam0 and then
%find a new estimate for the AR coefficients and white noise standard
%deviation using all "non-missing" trajectory.

if ~isempty(indx) %if there are missing data points

    %initialize unknowns
    arParam = arParam0;
    
    %initialize variable to be compared with arTol
    %to determine whether solution has been reached.
    arParamDiff = 10*arTol;
    
    %locate data points to be used in fitting (all but the first arOrder points + missing points)
    fitSet = indxP(arOrder+1:end);
    fitLength = length(fitSet);
    
    while arParamDiff > arTol %loop until required tolerance is reached
        
        %predict missing points using proposed AR model
%         [trajP,errFlag] = missPointARPred(traj,arParam);
        trajP = traj;
        for i = indx'
            trajP(i,1) = arParam*trajP(i-1:-1:i-arOrder,1);
        end
        
        %construct weighted matrix of previous points multiplying AR coefficients (on LHS of equation)
        %[size: fitLength by arOrder]
        lhsMat = zeros(fitLength,arOrder);
        for i=1:arOrder
            lhsMat(:,i) = trajP(fitSet-i,1)./trajP(fitSet,2);
        end
        
        %construct weighted vector of points used to determine parameters (on RHS of equation)
        %[size: fitLength by 1]
        rhsVec = trajP(fitSet,1)./trajP(fitSet,2);
        
        %get new estimate of arParam using least square fitting
        arParamNew = (lhsMat\rhsVec)';
        
        %compute variable to be compared with arTol
        arParamDiff = max(abs(arParamNew-arParam)); %maximum change in AR coefficients
        
        %update arParam
        arParam = arParamNew;
        
    end %(while arParamDiff > arTol)
    
    %check for causality of estimated model
    r = abs(roots([-arParam(end:-1:1) 1]));
    if ~isempty(find(r<=1.00001))
        disp('--arlsestimGaps: Warning: Predicted model (arParam) not causal!');
    end
    
    %get vector of weighted residuals
    epsilon = lhsMat*arParam' - rhsVec;
    
    %compute variance-covariance matrix
    varCovMat = inv(lhsMat'*lhsMat)*(epsilon'*epsilon/(fitLength-arOrder));
    
    %calculate standard deviation of white noise (note that nonweighted 
    %residuals must be used here)
    noiseSigma = std(epsilon.*trajP(fitSet,2));
    
    %predict missing points using estimated AR model
    [trajP,errFlag] = missPointARPred(traj,arParam);
        
else %if there are no missing data points

    arParam = [];
    varCovMat = [];
    noiseSigma = [];
    
end %(if ~isempty(indx))
