function [mu sigma] = getLifetimeDistributionMean(lftRes, histName)

if nargin<2
    histName = 'lftHist_A';
end

M = lftRes.(histName);
n = size(M,1);

t = lftRes.t;
dt = t(2)-t(1);
iMean = zeros(1,n);
for i = 1:n
    iM = M(setdiff(1:n,i),:);
%     iMean(i) = sum(mean(iM,1).*t*dt);
    iMean(i) = sum(M(i,:).*t*dt) / sum(M(i,:)*dt);
end
mu = mean(iMean);
sigma = std(iMean)/sqrt(n);
% sigma = std(iMean);
