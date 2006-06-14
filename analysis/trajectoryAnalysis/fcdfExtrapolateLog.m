function [logProb, isExtrapolated] = fcdfExtrapolateLog(fRatio, n1, n2)
%FCDFEXTRAPOLATELOG returns -log10(p) even for very small probabilities
%
% SYNOPSIS: [logProb, isExtrapolated] = fcdfExtrapolateLog(fRatio, n1, n2)
%
% INPUT fRatio: test statistic of the F-test (can be array)
%		n1/n2: degrees of freedom
%
% OUTPUT logProb: -log10(p), where p is the probability that the fRatio is 1. If p approaches 0, the logarithm will be extrapolated.
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
    % extrapolate the log probability. Calculate p-values in steps of 0.01
    % till we reach 1 to make sure that we get enough points where p is not
    % 0. Then, extrapolate log(-log(p))
    currentRatio = fRatio(i);
    currentN1 = n1(i);
    currentN2 = n2(i);
    
    ratioList = (currentRatio:0.01:1)';
    pList = fcdf(ratioList,currentN1, currentN2);
    % take loglog
    
    logLogList = log10(-log10(pList));
    logLogList(~isfinite(logLogList)) = 0;
    
    % find first 5 nonzero elements
    firstIdx = find(logLogList,5,'first');
    
    % linear fit
    linearFit = [ones(5,1),ratioList(firstIdx)]\logLogList(firstIdx);
    
    % extrapolate: y = lf(1) + lf(2)*x, therefore logProb =
    % 10^lf(1)*10^(lf(2)*x) = 10^(lf(1)+lf(2)*x)
    logProb(i) = 10^(linearFit(1) + linearFit(2) * currentRatio);
    
end

% reset warnings
warning(w);
