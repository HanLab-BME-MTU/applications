function [estimates, error] = ptFitCurve(xdata, ydata, tDrug)
% ptFitCurve calculates the best estimate for the coefficients a1, a2, t0 
% and t1 of the function:
%      t>t0 : a1
%      t0<=t<t1 : ((a2-a1)/(t1-t0))(t-t0)+a1 
%      t>=t1 : a2
%
% SYNOPSIS       ptFitCurve (xdata, ydata)
%
% INPUT          xdata : x-axis values
%                ydata : y-axis values
%                tDrug : timepoint where growth factor is put in the mix
%                
% OUTPUT         estimates : the best estimates for the coeff [a1 a2 t0 t1] 
%                error : is 1 if a curve cannot be fitted, 0 otherwise
%
% DEPENDENCIES   ptFitCurve uses { lsqcurvefit }
%                                  
%                ptFitCurve is used by { ptPlotEstimate }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Oct 04          Initial version

% Initialize error and estimates
error = 0;
estimates = [];

% Initialize the initial estimate matrix
guess = zeros(4,1);

% Provide initial estimate for b1
guess(1) = mean(ydata(1:tDrug));

% Calculate the standard deviation of b1
sigma = std(ydata(1:tDrug));

% Find the first value of ydata that is bigger than 3*sigma: this is
% going to be t0; we also want tDrug as a minimum
ind = find(ydata > guess(1)+(3*sigma));

% If index is empty it means that the graph is too flat for point > 3*sigma
% so let's try something smaller
if isempty(ind)    
   ind = find(ydata > guess(1)+(2*sigma));
   
   % If it's still empty, we give up because there won't be much to fit in
   % this data
   if isempty(ind) 
       error = 1;
       return;
   end
end

% Adjust guess if < tDrug
if ind(1) > tDrug
   guess(3) = ind(1);
else
   guess(3) = tDrug;
end

% Estimate the plateau at the end
guess(2) = mean(ydata(length(ydata)-tDrug:end));

% Estimate t1 by using the standard dev of the second plateau
sigma = std(ydata(length(ydata)-tDrug:end));
ind = find(ydata < guess(2)-(2*sigma));

% If index is empty it means that the graph is too flat for point < 2*sigma
% so let's try something smaller
if isempty(ind)
   ind = find(ydata < guess(2)-(1*sigma));
   
   % If it's still empty, we give up because there won't be much to fit in
   % this data
   if isempty(ind) 
       error = 1;
       return;
   end
end

% Adjust the t1 guess if < length(ydata)-tDrug
if ind(end) < length(ydata)-tDrug
   guess(4) = ind(end);
else
   guess(4) = length(ydata)-tDrug;
end

% Call lsqcurvefit with limits for t0 to come up with the best fit. The
% output is a vector with values [a1 a2 t0 t1]
estimates = lsqcurvefit(@expfun, guess, xdata, ydata, [-Inf; -Inf; tDrug; -Inf], ...
                        [Inf; Inf; length(ydata)-tDrug; Inf]);

% % expfun accepts curve parameters as inputs, and outputs a vector f. This
% function can be use by lsqcurvefit
function f = expfun(guess, x)
    % Check t0
    t0 = int16(guess(3));
    if t0 <= 0 
        t0 = 1;
        guess(3) = 1.0;
    end
    
    % Check t1
    t1 = int16(guess(4));
    if t1 > length(x) 
        t1 = length(x);
        guess(4) = double(length(x));
    end
    
    % Also make sure t1 does not come before t0
    if t1 < t0
        t1 = t0 + 1;
        guess(4) = guess(3) + 1.0;
    end
    
    % Retrieve vectors that we need
    xLeft = x(1:t0-1);
    xMiddle = x(t0:t1-1);
    xRight = x(t1:end);

    % Calculate the function value
    f1 = ones(1,length(xLeft)) .* guess(1);
    f2 = ((guess(2)-guess(1)) / (guess(4)-guess(3))) * (xMiddle-guess(3)) + guess(1);
    f3 = ones(1,length(xRight)) .* guess(2);
    f = [f1 f2 f3];
end

end  % function expfun