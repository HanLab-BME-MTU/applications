function [lhsVec,errFlag] = solveIterRefn(lhsMat,rhsVec,tol)
%solveIterRefn solves for x in Ax=b (matrix A, vector b) using iterative refinement
%
%SYNOPSIS [lhsVec,errFlag] = solveIterRefn(lhsMat,rhsVec,tol)
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
[nRow1,dummy] = size(lhsMat);
if dummy ~= nRow1
    disp('--solveIterRefn: "lhsMat" should be a square matrix!');
    errFlag = 1;
    lhsVec = [];
    return
end
[nRow2,dummy] = size(rhsVec);
if dummy ~= 1 || nRow2 ~= nRow1
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
    
%get initial estimate of solution
lhsVec = lhsMat\rhsVec;

%initialize variable to be compared with "tol"
change = 10*tol;

%iterate until specified tolerance is reached
while change > tol

    resid = rhsVec - lhsMat*lhsVec; %get residuals
    correction = lhsMat\resid; %solve for corrections
    lhsVec = lhsVec + correction; %update solution
    
    change = norm(correction)/norm(lhsVec); %get relative change to compare it with "tol"
    
end
