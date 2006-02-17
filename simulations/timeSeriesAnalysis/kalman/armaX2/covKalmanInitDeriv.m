function [stateCovMatDeriv00,errFlag] = covKalmanInitDeriv(F,FDeriv,...
    G,GDeriv,wnVariance,xOrder,covMat00,tolerance)
%COVKALMANINITDERIV calculates the partial derivatives of the initial covariance matrix used in Kalman recursions.
%
%SYNOPSIS [stateCovMatDeriv00,errFlag] = covKalmanInitDeriv(F,FDeriv,...
%    G,GDeriv,wnVariance,xOrder,covMat00,tolerance)
%
%INPUT  F         : Transition matrix in equation of state.
%       FDeriv    : Partial derivatives of transition matrix.
%       G         : Vector of influence of white noise process in
%                   equation of state.
%       GDeriv    : Vector of partial derivatives of G.
%       wnVariance: White noise variance.
%       xOrder    : Order of X component.
%       covMat00  : Initial state covariance matrix.
%       tolerance : Maximum acceptable difference between two
%                   consecutive estimates of the covariance matrix.
%
%OUTPUT stateCovMatDeriv00: Covariance matrix of initial state.
%       errFlag           : 0 if function executes normally, 1 otherwise.
%
%REMARKS This code solves for dP the equation 
%        dP = dF*P*F' + F*dP*F' + F*P*dF' + sigma^2*dG*G' + sigma^2*G*dG',
%        where dP is the unknown derivative of the covariance matrix, 
%        P is the known covariance matrix, F is the transition matrix, 
%        dF its derivative, G is the white noise influence, dG is its 
%        derivative and sigma^2 is the white noise variance.
%
%Khuloud Jaqaman, February 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stateCovMatDeriv00 = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if correct number of arguments was used when function was called
if nargin < nargin('covKalmanInitDeriv')
    disp('--covKalmanInitDeriv: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate derivative of covariance matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get size of covariance matrix
covSize = length(G);

%get size of matrix of derivatives
matSize = length(GDeriv);

%get number of parameters
numParam = matSize/covSize;

%set initial guess of covariance matrix to zero
stateCovMatDeriv00 = zeros(matSize,covSize);

for iparam = 1:numParam-(xOrder+1)

    %get relevant derivatives
    FDerivP = FDeriv((iparam-1)*covSize+1:iparam*covSize,1:end);
    GDerivP = GDeriv((iparam-1)*covSize+1:iparam*covSize);
    
    %set initial guess of derivative of covariance matrix to zero
    covDeriv = zeros(covSize);
    
    %initialize iteration variable
    doAgain = 1;

    %iterate until the derivative of the covariance matrix does not change
    %significantly from one iteration to another
    while doAgain

        %calculate new covariance matrix
        covDerivNew = FDerivP*covMat00*F' + F*covDeriv*F' + ...
            F*covMat00*FDerivP' + wnVariance*(GDerivP*G' + G*GDerivP');

        %calculate the relative difference between the old and new estimates of
        %the derivative
        relDiff = abs((covDerivNew-covDeriv)./covDerivNew);

        %indicate that iteration should stop if difference is smaller than
        %tolerance
        if max(relDiff) < tolerance
            doAgain = 0;
        end

        %assign new estimate to stateCovMat00
        covDeriv = covDerivNew;

    end %(while doAgain)

    %place derivative in proper place in big matrix
    stateCovMatDeriv00((iparam-1)*covSize+1:iparam*covSize,1:end) = covDeriv;
    
end %(for iparam = 1:numParam)


%%%%% ~~ the end ~~ %%%%%
