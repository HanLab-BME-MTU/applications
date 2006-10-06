function [logProb, isExtrapolated] = chi2cdfExtrapolateLog(chi2, n)
%CHI2CDFEXTRAPOLATELOG returns -log10(1-p) even for very small probabilities
%
% SYNOPSIS: [logProb, isExtrapolated] = chi2cdfExtrapolateLog(chi2, n)
%
% INPUT chi: test statistic of the chi2-test (can be array)
%		n: # of degrees of freedom.
%
% OUTPUT logProb: -log10(1-p), where 1-p is the probability that chi2 is n. If p approaches 0, the logarithm will be extrapolated.
%			isExtrapolated: true if the value has been extrapolated
%
% REMARKS
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: Jonas
% DATE: 12-Jun-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%================
%% TEST INPUT
%================

% check for input
if nargin ~= 2 || isempty(chi2) || isempty(n)
    error('chi2cdfExtrapolateLog needs three nonempty input arguments!')
end

% make inputs the correct size
[errorcode, chi2,n] = distchck(2,chi2,n);

if errorcode > 0
    error('chi2cdfExtrapolateLog requires non-scalar arguments to match in size.');
end


% preassign output, variables
[logProb,x12,x13] = deal(zeros(size(chi2)));

%================



%===================
%% GET PROBABILITY
%===================

% at sufficiently large chi2, or when 1-p drops below 1e-10, the logarithm
% becomes very linear. However, towards the end (1-p < 1e-14), chi2cdf returns
% "staircases". Thus, don't wait for the code to break down, but instead
% directly go for linear extrapolation if 1-p becomes less than 1e-10. To
% avoid unnecessary calls to chi2inv, we use an estimate for whether 1-p
% will be small.
% If 1-p is small, use [chi2inv(1e-12,n),-12] and [chi2inv(1e-13,n),-13] as
% two points to calculate a straight line. Since chi2inv breaks down above
% n=25, use estimates as the x-values, and calculate the
% log10(1-chi2cdf(x,n) for the two points.
% for n=1:15, use a cubic approximation for x, otherwise a linear one

% check for small p - estimate the chi2-value corresponding to 1-p=1e-12.
% since we're at it, do the same for 1e-13, too; it doesn't really increase
% computational time
bigNidx = n>10;

% use linear approximation for large n
x12(bigNidx) = 1.993*n(bigNidx) + 60.73;
x13(bigNidx) = 2.032*n(bigNidx) + 65.5;

% use cubic approximation for small n
x12(~bigNidx) = 0.0016*n(~bigNidx).^3-0.0939*n(~bigNidx).^2+3.8676*n(~bigNidx)+47.7630;
x13(~bigNidx) = 0.0016*n(~bigNidx).^3-0.0956*n(~bigNidx).^2+3.9467*n(~bigNidx)+52.2160;


% find big chi
bigChiIdx = chi2 > x12;

if any(bigChiIdx(:))
    % calculate y for the chis, and do linear extrapolation
    y12 = -log10(1-chi2cdf(x12(bigChiIdx),n(bigChiIdx)));
    y13 = -log10(1-chi2cdf(x13(bigChiIdx),n(bigChiIdx)));
    logProb(bigChiIdx) = y13 + ...
        (y13-y12)/(x13(bigChiIdx)-x12(bigChiIdx)) * ...
        (chi2(bigChiIdx)-x13(bigChiIdx));
end

if any(~bigChiIdx(:))
    % calculate logProb directly
    logProb(~bigChiIdx) = -log10(1-chi2cdf(chi2(~bigChiIdx),n(~bigChiIdx)));
end

% assign isExtrapolated - it's nothing else but our bigChiIdx
isExtrapolated = bigChiIdx;
