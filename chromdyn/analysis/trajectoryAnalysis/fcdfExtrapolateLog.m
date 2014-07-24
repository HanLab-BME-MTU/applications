function [logProb, isExtrapolated] = fcdfExtrapolateLog(fRatio, n1, n2)
%FCDFEXTRAPOLATELOG returns -log10(p) even for very small probabilities
%
% SYNOPSIS: [logProb, isExtrapolated] = fcdfExtrapolateLog(fRatio, n1, n2)
%
% INPUT fRatio: test statistic of the F-test (can be array)
%		n1/n2: degrees of freedom
%
% OUTPUT logProb: -log10(p), where p is the probability that the fRatio is
%                 1 given the degrees of freedom and the actual observed
%                 fRation for fRations <=1 (the code will take the inverse
%                 of fRatio if fRatio is above 1) If p approaches 0, the
%                 logarithm will be extrapolated.    
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
if nargin ~= 3 || isempty(fRatio) || isempty(n1) || isempty(n2)
    error('fcdfExtrapolateLog needs three nonempty input arguments!')
end

% make inputs the correct size
[errorcode, fRatio,n1,n2] = distchck(3,fRatio,n1,n2);

if errorcode > 0
    error('fcdfExtrapolateLog requires non-scalar arguments to match in size.');
end

% make sure ratio is smaller than 1, switch degrees of freedom if necessary
biggerIdx = fRatio > 1;
if any(biggerIdx)
    fRatio(biggerIdx) = 1./fRatio(biggerIdx);
    tmp = n1;
    n1(biggerIdx) = n2(biggerIdx);
    n2(biggerIdx) = tmp(biggerIdx);
end

% preassign output
logProb = zeros(size(fRatio));


% turn off log0-warning
w=warning;
warning('off','MATLAB:log:logOfZero');

%================



%===================
%% GET PROBABILITY
%===================

% try to get probability directly
fProb = fcdf(fRatio, n1, n2);

% find what we have to extrapolate
isExtrapolated = fProb == 0;


% set good probs
logProb(~isExtrapolated) = -log10(fProb(~isExtrapolated));

% for the loop, the idx has to be a row vector. Therefore, calculate it
% outside the loop already
extrapolatedIdx = find(isExtrapolated);

for i=extrapolatedIdx(:)'
    % extrapolate the log probability. For ratios below one, fcdf breaks
    % down at p=realmin, however, the function isn't linear as the ratio
    % approaches zero
    % x = 0.0001:0.0001:1;
    % figure,hold on, 
    % for j=1:100:10000,for i=j,
    % plot((x),-log10(fcdf(x,i,j)),...
    % 'Color',extendedColors( floor((j-1)/100)+1)),end,end
    % 
    % At ratios above one, log(-p) vs log(x) becomes linear, which we can
    % use for extrapolation.
    % Since fcdf already breaks down at p=1e-16 (becoming unstable before),
    % use the relation fcdf(x,i,j) = 1-fcdf(1/x,j,i) to calculate the
    % p-values.
    %
    % Unfortunately, there does not seem to be a simple way  to calculate 
    % where fcdf breaks down (or the level when e.g. p=1e-300):
    % for i=1:100
    % for j = 1:100
    % plist = -log10(fcdf(x,10^(0.05*i),10^(0.05*j)));
    % find the first entry that is NaN
    % logP(i,j) = x(find(isfinite(plist),1));
    % end
    % end
    % figure,surf(10.^(0.05:0.05:5),10.^(0.05:0.05:5),logP)
    %
    % Therefore, whenever NaN is encountered, the code will calculate
    % individually where the code breaks down, using the last two good
    % measurements (p1, p2) as a basis for the extrapolation
    % Better: Use the second and third good measurement, because the last
    % one can be iffy.
    %
    % -log10(pExtrapolated) = -log10(p2) + ...
    % (-log10(p2)+log10(p1))/(log10(x2)-log10(x1)) * ...
    % (log10(1/fRatio)-log10(x2))
    %
    % This will slow down the code somewhat, but it's still more stable
    % than trying finv, for example.

    % If we encounter a NaN in the calculation, we increase x in steps of
    % 0.001
    
    currentRatio = fRatio(i);
    currentN1 = n1(i);
    currentN2 = n2(i);
    
    % ratioList has to be the more finely spaced the closer to zero
    % currentRatio is.
    ratioList = [1:-0.0001:0.0002,1./(10000:1/currentRatio)]';
    pListLog = -log10(fcdf(ratioList,currentN1, currentN2));
    
    % find the first 5 non-Nan elements
    pIdx = find(isfinite(pListLog),5,'last');
    
    % read associated probabilities and ratios. Careful with the indices,
    % because low ratioList-values correspond to high x-values!
    lp1 = pListLog(pIdx(2));
    lp2 = pListLog(pIdx(1));
    x1 = 1/ratioList(pIdx(2));
    x2 = 1/ratioList(pIdx(1));
    
    % calculate logProb
    logProb(i) = lp2 + ...
    (lp2 - lp1)/(log10(x2)-log10(x1)) * ...
    (log10(1/currentRatio)-log10(x2));
    
end

% reset warnings
warning(w);
