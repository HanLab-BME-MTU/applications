function [estimates, sigma0, error] = ptFitLine(xdata, ydata)
% ptFitLine calculates the best estimate for the coefficient a of the function:
%      a * xdata = ydata
%
% SYNOPSIS       ptFitLine (xdata, ydata)
%
% INPUT          xdata : x-axis values
%                ydata : y-axis values
%                
% OUTPUT         estimates : the best estimates for the coeff [a1 a2 a3 a4] 
%                error : is 1 if a curve cannot be fitted, 0 otherwise
%
% DEPENDENCIES   ptFitLine uses {  }
%                                  
%                ptFitLine is used by { ptEstimateLine }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Oct 04          Initial version

% Initialize error and estimates
error = 0;
estimates = [];
guess = zeros(2,1);

% Initialize the initial estimate matrix
guess(1) = ydata(1);
guess(2) = 1;

% Call lsqcurvefit to come up with the best fit.
%[estimates,resnorm,residual,exitflag,output] = lsqcurvefit(@expfun, guess, xdata, ydata);

% Call leastMedianSquares to come up with the best fit, not taking into
% account any outliers
if size(xdata,2) == 1
    strData.x = xdata;
else
    strData.x = xdata';
end
if size(ydata,2) == 1
    strData.y = ydata;
else
    strData.y = ydata';
end
u0 = guess;
[estimates,goodRows,sigma0] = leastMedianSquares('y-(u(1)+(u(2)*x))', u0, [], strData);
                             


% expfun accepts curve parameters as inputs, and outputs a vector f. This
% function is used by lsqcurvefit
function f = expfun(guess, x)
    f = guess(1) + guess(2) * x;
end  % function expfun

end