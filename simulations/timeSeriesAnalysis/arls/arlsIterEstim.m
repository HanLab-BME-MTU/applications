function [arParam,varCovMat,residuals,noiseSigma,trajP,errFlag] = arlsIterEstim(traj,arOrder,arTol,numFuture)
%ARLSITERESTIM fits an AR model to a trajectory, and estimates missing data points, in an iterative fashion.
%
%SYNOPSIS [arParam,varCovMat,residuals,noiseSigma,trajP,errFlag] = arlsIterEstim(traj,arOrder,arTol,numFuture)
%
%INPUT  traj     : Trajectory to be modeled (with measurement uncertainty).
%                  Missing points should be indicated with NaN.
%       arOrder  : Order of proposed AR model.
%       arTol    : Change in AR coefficients below which it is safe to stop.
%       numFuture: Number of future time points to be used in prediction of
%                  a missing point (used in missPointARPred).
%
%OUTPUT arParam     : Parameters estimated by filling in the gaps and then using
%                     all points except for the missing ones.
%       varCovMat   : Variance-covariance matrix of arParam.
%       residuals   : Difference between measurements and model predictions.
%       noiseSigma  : Estimated standard deviation of white noise using AR
%                     model with parameters arParam and "non-missing" points.
%       trajP       : Trajectory with estimates of missing points and their uncertainties.
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, March 2004

errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= nargin('arlsIterEstim')
    disp('--arlsIterEstim: Incorrect number of input arguments!');
    errFlag  = 1;
    arParam = [];
    varCovMat = [];
    residuals = [];
    noiseSigma = [];
    trajP = []
    return
end

%check input data
if arOrder < 1
    disp('--arlsIterEstim: Variable "arOrder" should be >= 1!');
    errFlag = 1;
end
[trajLength,nCol] = size(traj);
if trajLength < arOrder
    disp('--arlsIterEstim: Length of trajectory should be larger than model order!');
    errFlag = 1;
end
if nCol ~= 2
    if nCol == 1
        traj = [traj ones(trajLength,1)];
    else
        disp('--arlsIterEstim: "traj" should have either 1 column for measurements, or 2 columns: 1 for measurements and 1 for measurement uncertainties!');
        errFlag = 1;
    end
end
if arTol <= 0
    disp('--arlsIterEstim: Variable "arTol" should be positive!');
    errFlag = 1;
end
if numFuture < 0 || numFuture > arOrder
    disp('--arlsIterEstim: numFuture should be between 0 and arOrder!')
    errFlag = 1;
end

%exit if there are problems in input data
if errFlag
    disp('--arlsIterEstim: Please fix input data!');
    arParam = [];
    varCovMat = [];
    residuals = [];
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

%check indx - algorithm cannot proceed if there are missing points for time <= arOrder.
if ~isempty(indx)
    if indx(1) <= arOrder
        disp('--arlsIterEstim: There are missing points for time points <= arOrder!');
        disp('  Please fix input data such that at least the first arOrder time points are available!');
        errFlag = 1;
        arParam = [];
        varCovMat = [];
        residuals = [];
        noiseSigma = [];
        trajP = [];
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate arParam0 using available data points that have enough previous points to
%predict them using the proposed order of the AR model. Note that if there
%are no missing data points, this is basically the solution!

[arParam0,varCovMat0,residuals0,noiseSigma0,errFlag] = arlsestim0(traj,arOrder);
if errFlag
    disp('--arlsIterEstim: Function arlsestim0 could not get initial estimate of parameters!');
    arParam = [];
    varCovMat = [];
    residuals = [];
    noiseSigma = [];
    trajP = [];
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if there are missing data points, predict them and then find a new
%estimate for the AR coefficients and white noise standard deviation
%using all available points. Repeat these two steps until coefficient
%estimates do not change significantly with further refinement.

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
        [trajP,errFlag] = missPointARPred(traj,arParam,numFuture);
        if errFlag
            disp('--arlsIterEstim: Function missPointARPred could not estimate missing points!');
            arParam = [];
            varCovMat = [];
            residuals = [];
            noiseSigma = [];
            trajP = [];
            return
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
        disp('--arlsIterEstim: Warning: Predicted model (arParam) not causal!');
    end
    
    %get vector of weighted residuals
    residuals = NaN*ones(trajLength,1);
    residuals(fitSet) = lhsMat*arParam' - rhsVec;
    
    %compute variance-covariance matrix
    varCovMat = inv(lhsMat'*lhsMat)*(residuals(fitSet)'*residuals(fitSet)/(fitLength-arOrder));
    
    %get nonweighted residuals
    residuals(fitSet) = residuals(fitSet).*traj(fitSet,2);
    
    %get standard deviation of white noise
    noiseSigma = std(residuals(fitSet));
    
    %predict missing points using estimated AR model
    [trajP,errFlag] = missPointARPred(traj,arParam,numFuture);
        
else %if there are no missing data points

    arParam = arParam0;
    varCovMat = varCovMat0;
    residuals = residuals0;
    noiseSigma = noiseSigma0;
    
end %(if ~isempty(indx))
