function [trajP,errFlag] = missPointTARPred(traj,vThresholds,delay,tarParam,future)
%MISSPOINTTARPRED estimates values of missing points in a trajectory assuming it is a TAR process
%
%SYNOPSIS [trajP,errFlag] = missPointTARPred(traj,vThresholds,delay,tarParam,future)
%
%INPUT  traj        : Trajectory to be modeled (with measurement uncertainty).
%                     Missing points should be indicated with NaN.
%       vThresholds : Column vector of thresholds, sorted in increasing order.
%       delay       : Time lag of value compared to vThresholds.
%       tarParam    : Coefficients of AR model.
%       future      : Logical variable with value '0' if no future points
%                     are used in estimation, and '1' if future points are used in
%                     estimation.
%
%OUTPUT trajP       : Trajectory with estimates of missing points and their uncertainties.
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, April 2004
%check if correct number of arguments were used when function was called

%initialize output
errFlag = 0;
trajP = [];

%check if correct number of arguments were used when function was called
if nargin ~= nargin('missPointTARPred')
    disp('--missPointTARPred: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check input data
if ~isempty(vThresholds)
    [nThresholds,nCol] = size(vThresholds);
    if nCol ~= 1
        disp('--missPointTARPred: "vThresholds" should be a column vector!');
    else
        if min(vThresholds(2:end)-vThresholds(1:end-1)) <= 0
            disp('--missPointTARPred: Entries in "vThresholds" should be sorted in increasing order, with no two elements alike!');
            errFlag = 1;
        end
    end
    if delay <= 0
        disp('--missPointTARPred: "delay" should be a positive integer!');
        errFlag = 1;
    end
else
    nThresholds = 0;
    delay = 1;
end
dummy = size(tarParam,1);
if dummy ~= nThresholds + 1
    disp('--missPointTARPred: Wrong number of rows in "tarParam"!');
    errFlag = 1;
else
    for i = 1:nThresholds+1
        tarOrder(i) = length(find(~isnan(tarParam(i,:))));
        r = abs(roots([-tarParam(i,tarOrder(i):-1:1) 1]));
        if ~isempty(find(r<=1.00001))
            disp('--missPointTARPred: Causality requires the autoregressive polynomial in each regime of the model not to have any zeros for z <= 1!');
            errFlag = 1;
        end
    end
end
[trajLength,nCol] = size(traj);
if trajLength < max(tarOrder)
    disp('--missPointTARPred: Length of trajectory should be larger than model order!');
    errFlag = 1;
end
if nCol ~= 2
    if nCol == 1
        traj = [traj ones(trajLength,1)];
    else
        disp('--missPointTARPred: "traj" should have either 1 column for measurements, or 2 columns: 1 for measurements and 1 for measurement uncertainties!');
        errFlag = 1;
    end
end
if future ~= 0 && future ~= 1
    disp('--missPointTARPred: "future" should be either 0 or 1!');
    errFlag = 1;
end

%exit if there are problems in input data
if errFlag
    disp('--missPointTARPred: Please fix input data!');
    return
end

%assign zeros to missing data points and ones to available measurements.
available = ~isnan(traj(:,1));

%get indices of missing points
indx = find(~available);

%check indx - algorithm cannot proceed if there are missing points for time <= max(max(tarOrder),delay).
if ~isempty(indx)
    if indx(1) <= max(max(tarOrder),delay)
        disp('--missPointTARPred: There are missing points for time points <= max(max(tarOrder),delay)!');
        disp('  Please fix input data such that all those time points are available!');
        errFlag = 1;
        return
    end
end

%replace NaN in entries of missing points with 0 to facilitate calculations later in the code
traj(indx,1) = 0;

%copy trajectory into trajP where missing points will be filled
trajP = traj;

%put +/- infinity at ends of thresholds vector
vThresholds = [-Inf; vThresholds; Inf];

%first predict all missing points by using the given TAR model and known previous points
for i = indx'
    level = find(((trajP(i-delay,1)>vThresholds(1:end-1)) + ... %determine level
        (trajP(i-delay,1)<=vThresholds(2:end))) == 2);
    trajP(i,1) = tarParam(level,1:tarOrder(level))*trajP(i-1:-1:i-tarOrder(level),1); %predict value
end
    
%if future = 1, estimate the missing points iteratively using values
%of points dependent on them as well.

%assign a level to each point
level = NaN*ones(trajLength,1);
for i = max(max(tarOrder),delay)+1:trajLength
    level(i) = find(((trajP(i-delay,1)>vThresholds(1:end-1)) + ...
        (trajP(i-delay,1)<=vThresholds(2:end))) == 2);
end

%copy trajP into trajP2
trajP2 = trajP;

%initialize variable indicating whether another iteration should be performed
doAgain = 1;
firstTime = 1;

if future
    
    while doAgain
        
        i0 = 1; %variable to loop over missing points
        
        while i0 <= length(indx) %go over all missing points
            
            i = i0;
            indxI = indx(i);
            missCluster = zeros(length(indx),1); %reserve memory for missCluster, delete missing points in previous cluster
            missCluster(i) = indxI; %first missing point in cluster
            
            if i0 ~= length(indx) %if missing point is not last missing point in trajectory
                
                numMissInt = 1;
                while numMissInt ~= 0
                    
                    indxI = indx(i); %time point at which data is missing
                    
                    jump = min(max(tarOrder),trajLength-indxI); %get farthest possible dependent point
                    while jump > tarOrder(level(indxI+jump)) && jump > 1 %determine dependence interval
                        jump = jump - 1;
                    end
                    
                    numMissInt = length(find(~available(indxI+1:indxI+jump))); %get number of missing points in this interval
                    missCluster(i+1:i+numMissInt) = indx(i+1:i+numMissInt); %add missing points to current cluster
                    
                    i = i + numMissInt; %update missing point index whose future points are checked
                    
                end
                
            end
            missCluster = missCluster(find(missCluster)); %remove all extra entries in array
            
            %number of missing points to be determined simultaneously
            numMiss = i-i0+1; 
            
            %number of equations to be used in determining missing points
            numEq = indx(i)-indx(i0)+1+max(tarOrder);
            
            %construct matrix on LHS multiplying unknowns
            %[size: numEq by numMiss]
            lhsMat = zeros(numEq,numMiss);
            for j = 1:numMiss %for each column
                start = missCluster(j)-missCluster(1)+1;
                finish = min(start+max(tarOrder),numEq);
                lhsMat(start,j) = 1;
                for j1 = start+1:finish
                    varNum = missCluster(j)+j1-start;
                    if varNum > trajLength
                        break
                    end
                    lhsMat(j1,j) = -tarParam(level(varNum),j1-start);
                end
            end
            
            %evaluate vector on RHS
            %[size: numEq by 1]
            rhsVec = zeros(numEq,1);
            for j=1:numEq
                varNum = indx(i0)+j-1;
                if varNum > trajLength
                    break
                end
                levelJ = level(varNum);
                rhsVecMult = [-1 tarParam(levelJ,1:tarOrder(levelJ))];
                rhsVec(j) = rhsVecMult*(available(varNum:-1:varNum-tarOrder(levelJ)).*...
                    trajP2(varNum:-1:varNum-tarOrder(levelJ),1));
            end
            
            %get rid of equations that do not involve any unknowns (such a situation
            %might arise due to different AR orders in different regimes
            lhsMat(find(isnan(lhsMat))) = 0; %first replace NaN with 0
            rowIndx = []; %then find rows in lhsMat with are all zero
            for j=1:numEq
                if ~isempty(find(lhsMat(j,:)))
                    rowIndx(end+1) = j; 
                end
            end
            lhsMat = lhsMat(rowIndx,:); %modify lhsMat
            rhsVec = rhsVec(rowIndx); %modify rhsVec
            NumEq = length(rowIndx); %update number of equations used in estimation
            
            %estimate missing data points
            trajP2(missCluster,1) = lhsMat\rhsVec;
            
            %get uncertainty in estimates. 
            %get vector of weighted residuals
            epsilon = lhsMat*trajP2(missCluster,1) - rhsVec;
            
            %compute variance-covariance matrix
            varCovMat = inv(lhsMat'*lhsMat)*(epsilon'*epsilon/(numEq-numMiss));
            
            %get uncertainty in estimates using varCovMat
            trajP2(missCluster,2) = sqrt(diag(varCovMat));
            
            %update variable indicating beginning of cluster
            i0 = i + 1;
            
        end %(while i0 <= length(indx))
        
        %assign a new level to each point
        levelP = NaN*ones(trajLength,1);
        for i = max(max(tarOrder),delay)+1:trajLength
            levelP(i) = find(((trajP2(i-delay,1)>vThresholds(1:end-1)) + ...
                (trajP2(i-delay,1)<=vThresholds(2:end))) == 2);
        end

        %compare new estimates of missing points to old ones, taking their
        %uncertainties into consideration
        if ~firstTime
            diff = abs(trajP2(indx,1)-trajP(indx,1)) - (trajP2(indx,2)+trajP(indx,2));
            if max(diff) < 0
                doAgain = 0;
            end
        end
        
        level = levelP;
        trajP = trajP2;
        firstTime = 0;
        
    end %(while levelChange)
    
end %(if future)
