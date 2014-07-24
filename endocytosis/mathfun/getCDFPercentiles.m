function p = getCDFPercentiles(t, cdfMat, pct)

N = size(cdfMat,1);
p = zeros(N, numel(pct));
for i = 1:N
    [~,idx] = unique(cdfMat(i,:));
    p(i,:) = interp1(cdfMat(i,idx), t(idx), pct);
end
