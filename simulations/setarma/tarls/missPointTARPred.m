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

errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= nargin('missPointTARPred')
    disp('--missPointTARPred: Incorrect number of input arguments!');
    errFlag  = 1;
    trajP = []
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
    trajP = [];
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
        trajP = [];
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

if future
    
    %assign a level to each point
    level = NaN*ones(trajLength,1);
    for i = max(max(tarOrder),delay)+1:trajLength
        level(i) = find(((trajP(i-delay,1)>vThresholds(1:end-1)) + ...
            (trajP(i-delay,1)<=vThresholds(2:end))) == 2);
    end

    %copy trajP into trajP2 where missing points will be estimated again
    trajP2 = trajP;
    
    i0 = 1; %variable to loop over all missing points
    
    while i0 <= length(indx) %go over all missing points
        
        i = i0;
        missCluster = zeros(length(indx),1); %delete missing points in previous cluster
        missCluster(i) = indx(i); %first missing point in cluster
        
        if i0 ~= length(indx)
            j = 1;
            while j ~= 0
                j = min(max(tarOrder),trajLength-i);
                while j > tarOrder(level(indx(i)+j)) && j > 1
                    j = j - 1;
                end
                while indx(i+1) <= indx(i)+j
                    i = i + 1;
                    missCluster(i) = indx(i);
                end
                if i == length(indx)
                    j = 0;
                end
            end
        end
        missCluster = missCluster(find(missCluster)); %remove all extra entries in array
        
        %number of missing points to be determined simultaneously
        numMiss = i-i0+1; 
        
        k = 1;
        for j = indx(i0):indx(i)+max(tarOrder)
            if max((j-indx(i0:i)>=0).*(j-indx(i0:i)<=tarOrder(level(j))))
                eqPart(k) = j;
                k = k + 1;
            end
        end
        
        %number of equation to be used in determining missing points
        numEq = length(eqPart);
        
%         %construct matrix on LHS multiplying unknowns
%         %[size: numEq by numMiss]
%         lhsMat = zeros(numEq,numMiss);
%         for j = 1:numMiss %for each column
%             start = indx(j+i0-1)-indx(i0)+1;
%             finish = min(start+arOrder,numEq);
%             lhsMat(start:finish,j) = [1 -arParam(1:finish-start)]';
%         end
%         
%         %evaluate vector on RHS
%         %[size: numEq by 1]
%         rhsVecMult = [-1 arParam];
%         rhsVec = zeros(numEq,1);
%         for j=1:numEq
%             varNum = indx(i0)+j-1;
%             if varNum > trajLength
%                 break
%             end
%             rhsVec(j) = rhsVecMult*(available(varNum:-1:varNum-arOrder).*...
%                 trajP(varNum:-1:varNum-arOrder,1));
%         end
%         
%         %in case end of trajectory is reached and less than the expected number
%         if j ~= numEq %of equations is available, 
%             lhsMat = lhsMat(1:j-1,:); %truncate lhsMat 
%             rhsVec = rhsVec(1:j-1); %and rhsVec accordingly
%         end
%         
%         %estimate missing data points
%         trajP(missCluster,1) = lhsMat\rhsVec;
%         
%         %get uncertainty in estimates. This is not possible if numFuture = 0
%         %(since all we have is one equation in one unknown)
%         if numFuture ~= 0
%             
%             %get vector of weighted residuals
%             epsilon = lhsMat*trajP(missCluster,1) - rhsVec;
%             
%             %compute variance-covariance matrix
%             varCovMat = inv(lhsMat'*lhsMat)*(epsilon'*epsilon/(numEq-numMiss));
%             
%             %get uncertainty in estimates using varCovMat
%             trajP(missCluster,2) = sqrt(diag(varCovMat));
%             
%         end
        
        %update variable indicating beginning of cluster
        i0 = i + 1;
        
        %update list of "available" points
        available(missCluster) = 1;
        
    end %(while i0 <= length(indx))
    
end %(if future)