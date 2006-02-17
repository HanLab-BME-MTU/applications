function [stateCovMat00,errFlag] = covKalmanInitGen(transitionMat,...
    procErrCov,wnVariance,tolerance)
%COVKALMANINITGEN calculates the covariance matrix of the initial state used in Kalman recursions.
%
%SYNOPSIS [stateCovMat00,errFlag] = covKalmanInitGen(transitionMat,...
%    procErrCov,wnVariance,tolerance)
%
%INPUT  transitionMat: Transition matrix in equation of state.
%       procErrCov   : Vector of influence of white noise process in
%                      equation of state (column vector).
%       wnVariance   : White noise variance.
%       tolerance    : Maximum acceptable difference between two
%                      consecutive estimates of the covariance matrix.
%
%OUTPUT stateCovMat00: Covariance matrix of initial state.
%       errFlag      : 0 if function executes normally, 1 otherwise.
%
%REMARKS This code solves the equation P = F*P*F' + sigma^2*G*G', where 
%        P is the unknown covariance matrix of the initial state, 
%        F is the transition matrix, G is the white noise influence
%        and sigma^2 is the white noise variance.
%
%Khuloud Jaqaman, February 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stateCovMat00 = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if correct number of arguments was used when function was called
if nargin < nargin('covKalmanInitGen')
    disp('--covKalmanInitGen: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculation of covariance matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get size of covariance matrix
matSize = length(procErrCov);

%set initial guess of covariance matrix to zero
stateCovMat00 = zeros(matSize);

%initialize iteration variable
doAgain = 1;

%iterate until the covariance matrix does not change significantly from one
%iteration to another
while doAgain

    %calculate new covariance matrix
    stateCovMatNew = transitionMat*stateCovMat00*transitionMat' + ...
        wnVariance*procErrCov*procErrCov';

    %calculate the relative difference between the old and new estimates of
    %the covariance matrix
    relDiff = abs((stateCovMatNew-stateCovMat00)./stateCovMatNew);
    
    %indicate that iteration should stop if difference is smaller than
    %tolerance
    if max(relDiff) < tolerance
        doAgain = 0;
    end

    %assign new estimate to stateCovMat00
    stateCovMat00 = stateCovMatNew;

end


%%%%% ~~ the end ~~ %%%%%
