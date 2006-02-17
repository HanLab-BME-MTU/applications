function [pacfV,errFlag] = pacf(gamma,maxLag)
%PACF calculates the partial autocorrelation function from the autocovariance function
%
%SYNOPSIS [pacfV,errFlag] = pacf(gamma)
%
%INPUT  gamma  : Autocovariance function of series or process, 
%                where gamma(i) is the autocovariance at lag i-1.
%       maxLag : Maximum lag at which PACF is calculated.
%
%OUTPUT pacfV  : Value of PACF upto maxLag.
%       errFlag: 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, February 2004

%initialize output
errFlag = 0;
pacfV = [];

%check if correct number of arguments were used when function was called
if nargin ~= nargin('pacf')
    disp('--pacf: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check input data
if maxLag <= 0
    disp('--pacf: Variable "maxLag" should be a positive integer!');
    errFLag = 1;
end
if length(gamma) < maxLag+1
    disp('--pacf: Length of "gamma" should be at least = "maxLag+1"!');
    errFlag = 1;
end
if errFlag
    disp('--pacf: Please fix input data!');
    return
end

%initialize pacfV
pacfV = zeros(maxLag+1,1);

%PACF at lag 1 always equals 1
pacfV(1) = 1;

for lag = 1:maxLag

    %get covariance matrix defined by Gamma(i,j) = gamma(i-j)
    covMat = diag(gamma(1)*ones(lag,1));
    for i = 1:lag-1
        entry = gamma(i+1)*ones(lag-i,1);
        covMat = covMat + diag(entry,i) + diag(entry,-i);
    end
    
    %solve equation covMat*soln=gamma
    soln = covMat\gamma(2:lag+1);
    
    %get partial autocorrelation function at current lag
    pacfV(lag+1) = soln(end);
    
end
