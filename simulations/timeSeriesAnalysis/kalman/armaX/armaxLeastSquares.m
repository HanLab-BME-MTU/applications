function [varCovMat,arParam,maParam,xParam,errFlag] = armaxLeastSquares(...
    trajOut,trajIn,wnVector,arOrder,maOrder,xOrder,constParam,wnVariance)
%ARMAXLEASTSQUARES estimates the ARMA coefficients and their variance-covariance matrix of a trajectory whose residuals are known.
%
%SYNOPSIS [varCovMat,arParam,maParam,xParam,errFlag] = armaxLeastSquares(...
%    trajOut,trajIn,wnVector,arOrder,maOrder,xOrder,constParam,wnVariance)
%
%INPUT  
%   Mandatory
%       trajOut   : Observations of output time series to be fitted. Either an
%                   array of structures trajOut(1:nTraj).observations, or a
%                   2D array representing one single trajectory.
%           .observations: 2D array of measurements and their uncertainties.
%                   Missing points should be indicated with NaN.
%       trajIn    : Observations of input time series. Either an
%                   array of structures trajIn(1:nTraj).observations, or a
%                   2D array representing one single trajectory.
%           .observations: 2D array of measurements and their uncertainties.
%                   Must not have any missing points.
%                   Enter as [] if there is no input series.
%       wnVector  : White noise series associated with a trajectory. Either
%                   an array of structures wnVector(1:nTraj).observations,
%                   or a 1D array of the WN series in 1 trajectory.
%           .observations: 1D array the white noise series in a trajectory.
%                   Enter as [] if there is no white noise series.
%       arOrder   : Order of AR part of process. Enter 0 if no AR part.
%   Optional
%       maOrder   : Order of MA part of process. Enter 0 if no MA part.
%                   Default: 0.
%       xOrder    : Order of X part of process. Enter -1 if no X part.
%                   Default: -1.
%       constParam: Set of constrained parameters. Constains 2 fields:
%           .ar      : 2D array. 1st column is AR parameter number and
%                      2nd column is parameter value. No need to input if
%                      there are no constraints on AR parameters.
%           .ma      : 2D array. 1st column is MA parameter number and
%                      2nd column is parameter value. No need to input if
%                      there are no constraints on MA parameters.
%           .x       : 2D array. 1st column is X parameter number and
%                      2nd column is parameter value. No need to input if
%                      there are no contraints on X parameters.
%                   Default: []
%       wnVariance: Estimated variance of white noise in process.
%                   Optional. Default: Zero.
%
%OUTPUT varCovMat : Variance-Covariance matrix of estimated coefficients:
%           .cofactorMat  : Cofactor matrix.
%           .posterioriVar: A posteriori estimate of residuals' variance.
%       arParam   : Estimated AR coefficients.
%       maParam   : Estimated MA coefficients.
%       xParam    : Estimated X coefficients.
%       errFlag   : 0 if function executes normally, 1 otherwise.
%
%REMARKS CONSTRAINED MINIMIZATION MUST BE UPDATED. DON'T USE!
%
%Khuloud Jaqaman, January 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varCovMat = [];
arParam = [];
maParam = [];
xParam = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if all mandatory input arguments were provided
if nargin < 4
    disp('--armaxLeastSquares: Incorrect number of input arguments!');
    errFlag = 1;
    return
end

%check trajOut and turn it into struct if necessary
if ~isstruct(trajOut)
    tmp = trajOut;
    clear trajOut
    trajOut.observations = tmp;
    clear tmp
elseif ~isfield(trajOut,'observations')
    disp('--armaxLeastSquares: Please input the trajOut in fields "observations"')
    errFlag = 1;
    return
end

%get number of trajectories supplied
numTraj = length(trajOut);

%check trajIn and turn it into struct if necessary
if isempty(trajIn) %if there is no input
    for i=1:numTraj
        trajIn(i).observations = [];
    end
else %if there is an input series
    if ~isstruct(trajIn) %turn into struct
        tmp = trajIn;
        clear trajIn
        trajIn.observations = tmp;
        clear tmp
    elseif ~isfield(trajIn,'observations')
        disp('--armaxLeastSquares: Please input trajIn in fields ''observations''!')
        errFlag = 1;
    end
end

%check wnVector and turn it into struct if necessary
if isempty(wnVector) %if there is no input
    for i=1:numTraj
        wnVector(i).observations = [];
    end
else %if there is an input series
    if ~isstruct(wnVector) %turn into struct
        tmp = wnVector;
        clear wnVector
        wnVector.observations = tmp;
        clear tmp
    elseif ~isfield(wnVector,'observations')
        disp('--armaxLeastSquares: Please input wnVector in fields ''observations''!')
        errFlag = 1;
    end
end

if arOrder < 0
    disp('--armaxLeastSquares: arOrder should be >= 0!');
    errFlag = 1;
end

if nargin < 5 || isempty(maOrder)
    maOrder = 0;
else
    if maOrder < 0
        disp('--armaxLeastSquares: maOrder should be >= 0!');
        errFlag = 1;
    elseif maOrder > 0 && isempty(wnVector)
        disp('--armaxLeastSquares: maOrder > 0 but there is no wnVector!');
        errFlag = 1;
    end
end

if nargin < 6 || isempty(xOrder)
    xOrder = -1;
else
    if xOrder < -1
        disp('--armaxLeastSquares: xOrder should be >= -1!');
        errFlag = 1;
    elseif xOrder > -1 && isempty(trajIn)
        disp('--armaxLeastSquares: xOrder > -1 but there is no trajIn!');
        errFlag = 1;
    end
end

%check parameter constraints
if nargin < 7 || isempty(constParam) %if no constraints were entered
    constParam = [];
else
    if isfield(constParam,'ar')
        nCol = size(constParam.ar,2);
        if nCol ~= 2
            disp('--armaxLeastSquares: constParam.ar should have 2 columns!');
            errFlag = 1;
        else
            if min(constParam.ar(:,1)) < 1 || max(constParam.ar(:,1)) > arOrder
                disp('--armaxLeastSquares: Wrong AR parameter numbers in constraint!');
                errFlag = 1;
            end
        end
    else
        constParam.ar = zeros(0,2);
    end
    if isfield(constParam,'ma')
        nCol = size(constParam.ma,2);
        if nCol ~= 2
            disp('--armaxLeastSquares: constParam.ma should have 2 columns!');
            errFlag = 1;
        else
            if min(constParam.ma(:,1)) < 1 || max(constParam.ma(:,1)) > maOrder
                disp('--armaxLeastSquares: Wrong MA parameter numbers in constraint!');
                errFlag = 1;
            end
        end
    else
        constParam.ma = zeros(0,2);
    end
    if isfield(constParam,'x')
        nCol = size(constParam.x,2);
        if nCol ~= 2
            disp('--armaxLeastSquares: constParam.x should have 2 columns!');
            errFlag = 1;
        else
            if min(constParam.x(:,1)) < 0 || max(constParam.x(:,1)) > maOrder
                disp('--armaxLeastSquares: Wrong X parameter numbers in constraint!');
                errFlag = 1;
            end
        end
    else
        constParam.x = zeros(0,2);
    end
end %(nargin < 7 || isempty(constParam) ... else ...)

%check white noise variance
if nargin < 8 || isempty(wnVariance) %if no white noise variance was entered
    wnVariance = 0;
else
    if wnVariance < 0
        disp('--armaxLeastSquares: White Noise Variance should be nonnegative!');
        errFlag = 1;
    end
end

%add/modify column for observational error of output
for i=1:numTraj
    [trajLength,nCol] = size(trajOut(i).observations);
    if nCol ~= 2
        if nCol == 1 %if no error is supplied, add a column of ones
            trajOut(i).observations = [trajOut(i).observations ...
                ones(trajLength,1)]; %assume that there is no observational error
        else % if there is more than 2 columns
            disp('--armaxLeastSquares: "trajOut.observations" should have either 1 column for measurements, or 2 columns: 1 for measurements and 1 for measurement uncertainties!');
            errFlag = 1;
        end
    end
end

for i=1:numTraj
    trajOut(i).observations(:,2) = sqrt(wnVariance + trajOut(i).observations(:,2).^2);
end

%exit if there are problems in input data
if errFlag
    disp('--armaxLeastSquares: Please fix input data!');
    return
end

%get maximum order and sum of orders
maxOrder = max([arOrder maOrder xOrder]);
sumOrder = arOrder + maOrder + xOrder + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameter estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fitLength = 0;
prevPoints = [];
observations = [];
epsilon = [];

%go over all trajectories
for i = 1:numTraj
    
    %obtain available points in this trajectory
    available = find(~isnan(trajOut(i).observations(:,1)));

    %find data points to be used in fitting
    fitSet = available(find((available(maxOrder+1:end)-available(1:end-maxOrder))...
        ==maxOrder) + maxOrder);
    fitLength1 = length(fitSet);
    fitLength = fitLength + fitLength1;
    
    %construct weighted matrix of previous output, white noise and input
    %multiplying  AR, MA and X coefficients, respectively
    %[size: fitLength1 by sumOrder]
    prevPoints1 = zeros(fitLength1,sumOrder);
    for j = 1:arOrder
        prevPoints1(:,j) = trajOut(i).observations(fitSet-j,1)...
            ./trajOut(i).observations(fitSet,2);
    end    
    for j = arOrder+1:arOrder+maOrder
        prevPoints1(:,j) = wnVector(i).observations(fitSet-j+arOrder,1)...
            ./trajOut(i).observations(fitSet,2);
    end
    for j = arOrder+maOrder+1:sumOrder
        prevPoints1(:,j) = trajIn(i).observations(fitSet-j+arOrder+maOrder+1,1)...
            ./trajOut(i).observations(fitSet,2);
    end
        
    %put points from all trajectories together in 1 matrix
    prevPoints = [prevPoints; prevPoints1];

    %form vector of weighted observations
    observations = [observations; trajOut(i).observations(fitSet,1)...
        ./trajOut(i).observations(fitSet,2)];

end %(for i = 1:numTraj)

if isempty(constParam) %if minimization in unconstrained

    %estimate ARMAX coefficients
    armaxCoef = (prevPoints\observations)';
    arParam = armaxCoef(1:arOrder);
    maParam = armaxCoef(arOrder+1:arOrder+maOrder);
    xParam  = armaxCoef(arOrder+maOrder+1:end);

    %get vector of weighted residuals
    epsilon = observations - prevPoints*armaxCoef';

    %calculate variance-covariance matrix
    varCovMat.cofactorMat = inv(prevPoints'*prevPoints);
    varCovMat.posterioriVar = epsilon'*epsilon/(fitLength-sumOrder);

else %if minimization is constrained

    %get cofactor matrix of unconstrained problem
    ucCofactMat = inv(prevPoints'*prevPoints);

    %get number of constraints
    numConst = length(constParam.ar(:,1)) + length(constParam.ma(:,1)) + ...
        length(constParam.x(:,1));
    constIndx = [constParam.ar(:,1); constParam.ma(:,1)+arOrder; ...
        constParam.x(:,1)+arOrder+maOrder+1];

    %construct matrix of constraints
    constMat = zeros(numConst,sumOrder);
    for i=1:numConst
        constMat(i,constIndx(i)) = 1;
    end

    %construct cofactor matrix of contrained problem
    dummy = inv(constMat*ucCofactMat*constMat');
    conCofactMat = [(ucCofactMat-ucCofactMat*constMat'*dummy*constMat* ...
        ucCofactMat) (ucCofactMat*constMat'*dummy); ...
        (dummy*constMat*ucCofactMat) (-dummy)];

    %combine observations vector with constraints
    observationsConst = [prevPoints'*observations; constParam.ar(:,2); ...
        constParam.ma(:,2); constParam.x(:,2)];

    %get solution
    lagrangeMult = (conCofactMat*observationsConst)';
    arParam = lagrangeMult(1:arOrder);
    maParam = lagrangeMult(arOrder+1:arOrder+maOrder);
    xParam = lagrangeMult(arOrder+maOrder+1:sumOrder);
    lagrangeMult = lagrangeMult(sumOrder+1:end);

    %get vector of weighted residuals
    armaxCoef = [arParam maParam xParam]';
    epsilon = observations - prevPoints*armaxCoef;

    %calculate variance-covariance matrix
    varCovMat.cofactorMat = conCofactMat(1:sumOrder,1:sumOrder);
    varCovMat.posterioriVar = epsilon'*epsilon/(fitLength-sumOrder);

end

%%%%% ~~ the end ~~ %%%%%
