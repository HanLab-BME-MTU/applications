function [lhsVec,errFlag] = lsIterRefn(lhsMat,rhsVec,tol)%solveIterRefn solves for x in Ax=b (matrix A, vector b) using iterative refinement
%LSITERREFN solves for x in the overdetermined problem Ax=b using iterative refinement
%
%SYNOPSIS [lhsVec,errFlag] = lsIterRefn(lhsMat,rhsVec,tol)
%
%INPUT  lhsMat : Matrix multiplying unknown on left hand side of equation.
%       rhsVec : Vector on right hand side of equation.
%       tol    : Tolerance at which calculation stops.
%
%OUTPUT lhsVec : Unknown vector on left hand side of equation, to be solved for.
%       errFlag: 0 if function executes normally, 1 otherwise.
%
%COMMENTS Right now, I'm using matrix left division each step. The code can
%         be written in a more efficient way emplying the same matrix 
%         triangularization all the time.
%
%Khuloud Jaqaman, March 2004

errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= nargin('solveIterRefn')
    disp('--solveIterRefn: Incorrect number of input arguments!');
    errFlag  = 1;
    arParam = [];
    noiseSigma = [];
    return
end

%check input data
[nRow,nCol] = size(lhsMat); %get # of data points, # of unknowns
if nCol >= nRow
    disp('--solveIterRefn: "lhsMat" should have more rows than columns!');
    errFlag = 1;
    lhsVec = [];
    return
end
[dummy1,dummy2] = size(rhsVec);
if dummy2 ~= 1 || dummy1 ~= nRow
    disp('--solveIterRefn: "rhsVec" should be a column vector of length equal to number of rows in "lhsMat"!');
    errFlag = 1;
end
if tol <= 0
    disp('--solveIterRefn: "tol" should be a small, positive number!');
    errFlag = 1;
end
%exit if there are problems in input data
if errFlag
    lhsVec = [];
    return
end

%initialize solution and residual vectors
lhsVec = zeros(nCol,1);
resid = zeros(nRow,1);

%initialize variable to be compared with "tol"
change = 10*tol;

%iterate until specified tolerance is reached
while change > tol
    
    %solve for corrections to solution and residuals
    tempLhsMat = [eye(nRow) lhsMat; lhsMat' zeros(nCol)]; %matrix on LHS
    tempRhsVec = [rhsVec; zeros(nCol,1)] - tempLhsMat*[resid; lhsVec]; %vector on RHS
    correction = tempLhsMat\tempRhsVec; %get correction vector
    
    %update solution and residuals
    resid = resid + correction(1:nRow);
    lhsVec = lhsVec + correction(nRow+1:end);
    
    %get relative change to compare it to "tol"
    change = norm(correction)/norm([resid; lhsVec]);
    
end
