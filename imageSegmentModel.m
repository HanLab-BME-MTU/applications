function Im = imageSegmentModel(imSize, sigmaPSF, segmentParams)

[n,p] = size(segmentParams);

if p ~= 5
    error('Invalid number of segment parameters.');
end

% xy-coordinates of the segment centers
xC = segmentParams(:, 1);
yC = segmentParams(:, 2);

% average amplitude
A = segmentParams(:, 3);

% length of the segment
L = segmentParams(:,4);
L2 = L / 2;

% Orientation of the segment
t = segmentParams(:,5);

% Hypothenuse length, corresponding to the half-length diagonal of a
% 2*(L2+d) long by 2d wide rectangle surrounding the segment.
d = 5 * sigmaPSF;
lh = sqrt(d.^2 + (L2 + d).^2);

% Angle between a rectangle border and a diagonal of the rectangle
at = atan(d ./ (L2 + d));

s1 = repmat([1 1 -1 -1], n, 1);
s2 = repmat([1 -1 1 -1], n, 1);

% xy-coordinates of the 4 rectangle's corners.
x = repmat(xC,1,4) + s1 .* cos(repmat(t,1,4) + s2 .* repmat(at,1,4)) .* repmat(lh,1,4);
y = repmat(yC,1,4) + s1 .* sin(repmat(t,1,4) + s2 .* repmat(at,1,4)) .* repmat(lh,1,4);

xMin = max(min(floor(x),[],2), ones(n,1));
xMax = min(max(ceil(x),[],2), ones(n,1) * imSize(2));
yMin = max(min(floor(y),[],2), ones(n,1));
yMax = min(max(ceil(y),[],2), ones(n,1) * imSize(1));

% Generate the image segments
Im = zeros(imSize);

for i = 1:n
    xRange = xMin(i):sign(xMax(i)-xMin(i)):xMax(i);
    yRange = yMin(i):sign(yMax(i)-yMin(i)):yMax(i);
    
    S = dLSegment2D(xRange, yRange, xC(i), yC(i), A(i),...
        sigmaPSF, L(i), t(i));

    Im(yRange,xRange) = Im(yRange,xRange) + S;
end