function [testStatistic,errFlag] = globalCoefTestStat(armaxParam1,armaxParam2,...
    varCovMatT1,varCovMatT2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

testStatistic = [];
errFlag =  0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check input data
if nargin ~= 4 
    disp('--globalCoefTestStat: Wrong number of input arguments!');
    errFlag = 1;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Test-statistic calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate vector of differences in coefficients
diffM = armaxParam1 - armaxParam2;

%calculate variance-covariance matrix of difference vector
diffV = varCovMatT1 + varCovMatT2;

%compute testStatistic
% testStatistic = diffM*(diffV\diffM')/length(diffM);
testStatistic = diffM*(diffV\diffM');


%%%%% ~~ the end ~~ %%%%%
