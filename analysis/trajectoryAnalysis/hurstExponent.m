function hurstParms = hurstExponent(trajectory, verbose)
%HURSTEXPONENT finds trajectory descriptors based on the correlation function
%
% The Hurst parameters are used in solid state physics to describe rough
% interfaces by fitting the autocorrelation function corr(Lag) with 
% sigma^2 * exp-{(Lag/xi)^2h}
%
% SYNOPSIS hurstParms = hurstExponent(trajectory)
%
% INPUT trajectory: any 1-D vector
%       verbose   : (optional) whether or not to plot resulting fit [{0}/1]
%
% OUTPUT hurstParms [sigma, xi, h] 
%
%c: jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%==============
% TEST INPUT
%==============

% defaults
defVerbose = 0;

% not much to do here, just make it a col-vector
trajectory = trajectory(:);

% check verbose
if nargin > 1 && ~isempty(verbose)
    % take input
else
    verbose = defVerbose;
end
    

%==============


%===============================
% CALCULATE FIRST ESTIMATE
%===============================

% calculate correlation, remove linear trends
lag = [0:20]';
lagX = [1:20]';
corr = autoCorr(trajectory, lag(end), 0);
corr = abs(corr);
corrX = abs(corr(2:end));

% calculate sigma
sigma = nanStd(trajectory);

% linear fit: y = a*x + b with
% y = ln(2*ln(sigma) - ln(Corr))
% a = 2*h
% x = lag
% b = -2*h*ln(xi)
% however, corr is already normalized, so no need for sigma

% put into form A*x = B with x = [b, a]

B = log( - log(corrX));
A = [ones(lagX(end),1), log(lagX)];

% linear fit
x1 = A\B;
%x2 = linearLeastMedianSquares(A,B);

% calculate Hurst parms
h1 = x1(2)/2;
xi1 = exp(-x1(1)/x1(2));
%h2 = x2(2)/2;
%xi2 = exp(-x2(1)/x2(2));

if verbose
fh=figure;
plot(lag, corr,'k',lag, exp(-(lag/xi1).^(2*h1)),'g');%,...
    %lag,exp(-(lag/xi2).^(2*h2)),'b')
end

%===============================


%========================
% FINAL FIT
%========================

% use lsqcurvefit because we can easily input ydata
options = optimset('Display','off');
hurst = lsqcurvefit(@fitfun, [xi1, h1], lag, corr, [],[], options)

if verbose
% plot
figure(fh);
hold on
plot(lag, exp(-(lag/hurst(1)).^(2*hurst(2))),'r');
end

% output
hurstParms = [sigma,xi1,h1;sigma,hurst];

%========================



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SUBFUNCTION FITFUN
function estimate = fitfun(hurst,lag)
% fitfun calculates the hurst approximation to the autocorrelation function

xi = hurst(1);
h  = hurst(2);

estimate = exp(-(lag/xi).^(2*h));
