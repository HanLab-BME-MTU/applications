function [estimates, error] = ptFitCurve(xdata, ydata, tDrug, slope)
% ptFitCurve calculates the best estimate for the coefficients a1, a2, t0 
% and t1 of the function:
%      a1 + a2 * tanh(a3*(x-a4))
%
% SYNOPSIS       ptFitCurve (xdata, ydata, tDrug, slope)
%
% INPUT          xdata : x-axis values
%                ydata : y-axis values
%                tDrug : timepoint where growth factor is put in the mix
%                slope : positive (>0) or negative (<0) slope expected
%                
% OUTPUT         estimates : the best estimates for the coeff [a1 a2 a3 a4] 
%                error : is 1 if a curve cannot be fitted, 0 otherwise
%
% DEPENDENCIES   ptFitCurve uses { leastMedianSquares }
%                                  
%                ptFitCurve is used by { ptPlotEstimate }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Oct 04          Initial version
% Andre Kerstens        Dec 04          Added slope parameter to be able to
%                                       handle neg and pos slopes (user can choose)
% Andre Kerstens        Mar 05          Changed default fit function to
%                                       tangens hyperbolicus (tanh) and fit
%                                       is made using leastMedianSquares

% Initialize error and estimates
error = 0;
estimates = [];

if ~exist('slope','var')
   slope = 1;
end

if ~exist('tDrug','var')
   tDrug = 35;
end

% Initialize the initial estimate matrix
guess = zeros(4,1);

% Plateaus
if slope >= 0
    platLow = mean(ydata(1:tDrug));
    platHigh = mean(ydata(length(ydata)-tDrug:end));
else
    platLow = mean(ydata(length(ydata)-tDrug:end));
    platHigh = mean(ydata(1:tDrug));
end

% Provide initial estimate for a2 
guess(2) = (platHigh - platLow) / 2;

% Provide initial estimate for a1 
guess(1) = guess(2) + platLow;
                                         
% Calculate the standard deviation of the lower plateau
if slope >= 0
    sigmaLow = std(ydata(1:tDrug));
else
    sigmaLow = std(ydata(length(ydata)-tDrug:end));
end

% Find the first value of ydata that is bigger/smaller than 3*sigma: this is
% going to be the estimate for t0; we also want tDrug as a minimum
if slope >= 0
   ind = find(ydata > platLow+(3*sigmaLow));
else
   ind = find(ydata < platLow-(3*sigmaLow));
end

% If index is empty it means that the graph is too flat for point >/< 3*sigma
% so let's try something smaller
if isempty(ind)  
   if slope >= 0
      ind = find(ydata > platLow+(2*sigmaLow));
   else
      ind = find(ydata < platLow-(2*sigmaLow));
   end
   
   % If it's still empty, we give up because there won't be much to fit in
   % this data
   if isempty(ind) 
       error = 1;
       return;
   end
end

% Adjust guess if < tDrug
if ind(1) > tDrug
   tDrug = ind(1);
end

% Estimate t1 by using the standard dev of the second plateau
if slope >= 0
    sigmaHigh = std(ydata(length(ydata)-tDrug:end));
else
    sigmaHigh = std(ydata(1:tDrug));
end
if slope >= 0
   ind = find(ydata < platHigh-(2*sigmaHigh));
else
   ind = find(ydata > platHigh+(2*sigmaHigh));
end

% If index is empty it means that the graph is too flat for point < 2*sigma
% so let's try something smaller
if isempty(ind)
   if slope >= 0
      ind = find(ydata < platHigh-(1*sigmaHigh));
   else
      ind = find(ydata > platHigh+(1*sigmaHigh));
   end
   
   % If it's still empty, we give up because there won't be much to fit in
   % this data
   if isempty(ind) 
       error = 1;
       return;
   end
end

% Adjust the t1 guess if < length(ydata)-tDrug
if ind(end) < length(ydata)-tDrug
   tDrugEnd = ind(end);
else
   tDrugEnd = length(ydata)-tDrug;
end

% Estimate for the slope value a3
guess(3) = (platHigh - platLow) / (tDrugEnd - tDrug);

% Estimate for a4 (shift from 0)
guess(4) = round(tDrug + (tDrugEnd - tDrug)/2);

% Call lsqcurvefit with limits for t0 to come up with the best fit. The
% output is a vector with values [a1 a2 a3 a4]
%estimates = lsqcurvefit(@expfun, guess, xdata, ydata);

% Call leastMedianSquares to come up with the best fit, not taking into
% account any outliers
strData.x = xdata';
strData.y = ydata';
u0 = guess;
options = optimset('MaxFunEvals',500);
[estimates,goodRows,sigma0] = leastMedianSquares('y-(u(1)+u(2)*tanh(u(3)*(x-u(4))))',...
                                         u0, options, strData);
                                     
% expfun accepts curve parameters as inputs, and outputs a vector f. This
% function is used by lsqcurvefit
function f = expfun(guess, x)
    f = guess(1) + guess(2)*tanh(guess(3)*(x - guess(4)));
end

end  % function expfun