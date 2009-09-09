function CC = createLabelsFromWindows(winPoly, N, M)

[nBands,nWindows] = size(winPoly);

CC = zeros(N, M);

for m = 1:nBands
    for n = 1:nWindows
        if ~isempty(winPoly(m,n).outerBorder) && ~isempty(winPoly(m,n).innerBorder) && ...
           ~max(isnan(winPoly(m,n).outerBorder(:))) && ~max(isnan(winPoly(m,n).innerBorder(:)))
            l = sub2ind([nBands, nWindows], m, n);
            
            X = [winPoly(m,n).outerBorder(1, :), winPoly(m,n).innerBorder(1, end:-1:1)];
            Y = [winPoly(m,n).outerBorder(2, :), winPoly(m,n).innerBorder(2, end:-1:1)];

            BW = poly2mask(X, Y, N, M);
            
            CC(BW == 1) = l;
        end
    end
end

end