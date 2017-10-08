function xcvec = nanXcorrelation(x, y, maxLag)
% nanXcorrelation Compute cross correlations between x and y vectors in a
% fashion that can handle many NaN's. It computes a cross product of
% centered x and y, and then normalize it based nanstd and effective sample
% sizes of x and y vectors. If there is no NaN value, then it coincides
% with xcorr.m function.
%
% Jungsik Noh, 2016/10


if numel(x) ~= numel(y)
    error('X and Y have different number of observations');
elseif numel(x) < maxLag
    error('maxLag should be less than numel of the input vector')
end

nx = numel(x);
xcvec = nan(1, 2*maxLag+1);

x = x(:);
y = y(:);

normFactor = nanstd(x)*sqrt(sum(isfinite(x))-1) * nanstd(y)*sqrt(sum(isfinite(y))-1);

if isfinite(normFactor)
    % Idea: cross prod with nansum
    x = x - nanmean(x);
    y = y - nanmean(y);

    x1 = nan(2*maxLag + nx, 1);
    x1(maxLag+1:maxLag+nx) = x;

    for i = 1:(2*maxLag+1)
        a = x1(i:i+nx-1);
        if all(isnan(a .* y))
            xcvec(i) = NaN;
        else
            xcvec(i) = sum(a .* y, 'omitnan');
        end
    end

    xcvec = xcvec./ normFactor;
end
    


end
