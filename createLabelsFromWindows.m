function CC = createLabelsFromWindows(winPoly, N, M)

[nBands, nWindows] = size(winPoly);

CC = zeros(N, M);
for n = 1:nWindows
    e = cellfun(@(x) ~(isempty(x) || any(isnan(x(:)))), {winPoly(:,n).outerBorder}) & ...
        cellfun(@(x) ~(isempty(x) || any(isnan(x(:)))), {winPoly(:,n).innerBorder});

    mStart = find(e, 1, 'first');
    mEnd = find(e, 1, 'last');
    
    % This is a little patch to cope non-contigous windows in a slice.
    % Should not happen anymore.
    if ~min(e(mStart:mEnd))
        continue;
    end
    
    X = winPoly(mStart, n).outerBorder(1, :);
    Y = winPoly(mStart, n).outerBorder(2, :);
        
    for m = mStart+1:mEnd-1
       X = horzcat(X, winPoly(m, n).outerBorder(1, end));
       Y = horzcat(Y, winPoly(m, n).outerBorder(2, end));
    end

    X = horzcat(X, winPoly(mEnd, n).innerBorder(1, end:-1:1));
    Y = horzcat(Y, winPoly(mEnd, n).innerBorder(2, end:-1:1));
        
    for m = mEnd-1:-1:mStart+1
       X = horzcat(X, winPoly(m, n).outerBorder(1, 1));
       Y = horzcat(Y, winPoly(m, n).outerBorder(2, 1));
    end
    
    BW =  poly2mask(X, Y, N, M);
    
    CC(BW == 1) = n;    
end

end