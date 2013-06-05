%[K] = corrMatFromCov(C) returns the correlation matrix given a covariance matrix

function [K] = corrMatFromCov(C)
n = size(C,1);
K = zeros(n,n);

idx = pcombs(1:n);
i = idx(:,1);
j = idx(:,2);
ij = i+n*(j-1);
ii = i+n*(i-1);
jj = j+n*(j-1);

K(ij) = C(ij) ./ sqrt(C(ii).*C(jj));

% remaining components are redundant
K = K + K';
K(sub2ind([n n], 1:n, 1:n)) = 1;
