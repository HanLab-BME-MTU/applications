function fittedCorVec = smoothingSplineCorMap(corMap)
% smoothingSplineCorMap Compute an averaged correlation curve over windows
% by applying smoothingspline via fit.m function. The smoothingParam is
% chosen to be twice of the default one. This function handles NaN
% correlations by simply removing in the fit.
% Jungsik Noh, 2017/03/31

xl = 1:size(corMap, 2);
yy = corMap(:);
tmpA = repmat(xl, size(corMap, 1), 1);
xx = tmpA(:);

% remove nan's
xx = xx(~isnan(yy));
yy = yy(~isnan(yy));

[~,~,tmpfit] = fit(xx, yy, 'smoothingspline');
smParSp = min(1, tmpfit.p*2);   % make the fit less smooth than the optimal

fitsp = fit(xx, yy, 'smoothingspline', 'SmoothingParam', smParSp);
fittedCorVec = fitsp(xl);

fittedCorVec = reshape(fittedCorVec, 1, []);

end
