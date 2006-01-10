function [trajP,errFlag] = missPointARPred(traj,arParam,numFuture)
%MISSPOINTARPRED estimates values of missing points in a trajectory assuming it is an AR(p) process
%
%SYNOPSIS [trajP,errFlag] = missPointARPred(traj,arParam,numFuture)
%
%INPUT  traj     : Trajectory to be modeled (with measurement uncertainty).
%                  Missing points should be indicated with NaN.
%       arParam  : Coefficients of AR model.
%       numFuture: Number of future time points to be used in prediction of
%                  a missing point.
%
%OUTPUT trajP       : Trajectory with estimates of missing points and their uncertainties.
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, March 2004
%check if correct number of arguments were used when function was called

errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= nargin('missPointARPred')
    disp('--missPointARPred: Incorrect number of input arguments!');
    errFlag  = 1;
    trajP = [];
    return
end

%get order of AR model
arOrder = length(arParam);

%check input data
[trajLength,nCol] = size(traj);
if trajLength < arOrder
    disp('--missPointARPred: Length of trajectory should be larger than model order!');
    errFlag = 1;
end
if nCol ~= 2
    if nCol == 1
        traj = [traj ones(trajLength,1)];
    else
        disp('--missPointARPred: "traj" should have either 1 column for measurements, or 2 columns: 1 for measurements and 1 for measurement uncertainties!');
        errFlag = 1;
    end
end
r = abs(roots([-arParam(end:-1:1) 1]));
if ~isempty(find(r<=1.00001))
    disp('--missPointARPred: Input model must be causal!');
    errFlag = 1;
end
if numFuture < 0 || numFuture > arOrder
    disp('--missPointARPred: numFuture should be between 0 and arOrder!');
    errFlag = 1;
end

%exit if there are problems in input data
if errFlag
    disp('--missPointARPred: Please fix input data!');
    trajP = [];
    return
end

%assign zeros to missing data points and ones to available measurements.
available = ~isnan(traj(:,1));

%get indices of missing points
indx = find(~available);

%check indx - algorithm cannot proceed if there are missing points for time <= arOrder.
if ~isempty(indx)
    if indx(1) <= arOrder
        disp('--missPointARPred: There are missing points for time points <= arOrder!');
        disp('  Please fix input data such that at least the first arOrder time points are available!');
        errFlag = 1;
        trajP = [];
        return
    end
end

%replace NaN in entries of missing points with 0 to facilitate calculations later in the code
traj(indx,1) = 0;

%copy trajectory into trajP where missing points will be filled
trajP = traj;

i0 = 1; %variable to loop over all missing points

while i0 <= length(indx) %go over all missing points
    
    i = i0;
    missCluster = zeros(length(indx),1); %delete missing points in previous cluster
    missCluster(i) = indx(i); %first missing point in cluster
    
    if i0 ~= length(indx) %if first missing point in cluster is not last missing point in trajectory,
        while indx(i+1)-indx(i) <= numFuture %check next missing point to see if it is less 
            i = i + 1;                       %than or equal to numFuture points away.
            missCluster(i) = indx(i); %in that case, put next point in same cluster
            if i == length(indx)
                break
            end
        end
    end %otherwise, this missing point is alone in the cluster
    missCluster = missCluster(find(missCluster)); %remove all extra entries in array
    
    %number of missing points to be determined simultaneously
    numMiss = i-i0+1; 
    
    %number of equation to be used in determining missing points
    numEq = indx(i)-indx(i0)+1+numFuture;
    
    %construct matrix on LHS multiplying unknowns
    %[size: numEq by numMiss]
    lhsMat = zeros(numEq,numMiss);
    for j = 1:numMiss %for each column
        start = indx(j+i0-1)-indx(i0)+1;
        finish = min(start+arOrder,numEq);
        lhsMat(start:finish,j) = [1 -arParam(1:finish-start)]';
    end
    
    %evaluate vector on RHS
    %[size: numEq by 1]
    eqBreak = 0;
    rhsVecMult = [-1 arParam];
    rhsVec = zeros(numEq,1);
    for j=1:numEq
        varNum = indx(i0)+j-1;
        if varNum > trajLength
            eqBreak = 1;
            break
        end
        rhsVec(j) = rhsVecMult*(available(varNum:-1:varNum-arOrder).*...
            trajP(varNum:-1:varNum-arOrder,1));
    end
    
    %in case end of trajectory is reached and less than the expected number
    if eqBreak %of equations is available, 
        lhsMat = lhsMat(1:j-1,:); %truncate lhsMat 
        rhsVec = rhsVec(1:j-1); %and rhsVec accordingly
    end
    
    %estimate missing data points
    trajP(missCluster,1) = lhsMat\rhsVec;
    
    %get uncertainty in estimates. This is not possible if numFuture = 0
    %(since all we have is one equation in one unknown)
    if numFuture ~= 0
        
        %get vector of weighted residuals
        epsilon = lhsMat*trajP(missCluster,1) - rhsVec;
        
        %compute variance-covariance matrix
        varCovMat = inv(lhsMat'*lhsMat)*(epsilon'*epsilon/(numEq-numMiss));
        
        %get uncertainty in estimates using varCovMat
        trajP(missCluster,2) = sqrt(diag(varCovMat));
    
    end
    
    %update variable indicating beginning of cluster
    i0 = i + 1;
    
    %update list of "available" points
    available(missCluster) = 1;
    
end
