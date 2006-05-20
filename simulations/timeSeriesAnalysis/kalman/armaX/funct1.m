function fc = funct1(n,xc)
%FUNCT1 is called by NAG's e04jaf to calculate -2ln(likelihood)
%
%SYNOPSIS fc = funct1(n,xc)
%
%INPUT  n : Number of parameters.
%       xc: Values of parameters.
%
%OUTPUT fc: -2ln(likelihood) evaluated at given set of parameters.
%
%Khuloud Jaqaman, September 2004

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fc = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin ~= 2
    disp('--funct1: Incorrect number of input arguments!');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Likelihood calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get additional input parameters
load funct1Input;

%evaluate -2ln(likelihood)
fc = neg2LnLikelihoodX(xc,prob);


%%%%% ~~ the end ~~ %%%%%

