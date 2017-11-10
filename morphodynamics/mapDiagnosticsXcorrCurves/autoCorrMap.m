function acfMat1 = autoCorrMap(map, varargin)
% autoCorrMap Compute auto-correlations over time for each window (over
% columns for each row) while handling NaN's.
% Jungsik Noh, 2016/10/29

tmax = size(map, 2);
wmax = size(map, 1);
maxLag0 = floor(tmax/4);

ip = inputParser;
ip.addRequired('map', @ismatrix)
ip.addParameter('maxLag', maxLag0);

ip.parse(map, varargin{:});
maxLag = ip.Results.maxLag;

acfMat1 = zeros(wmax, maxLag+1);

for w = 1:wmax
    acf0 = nanXcorrelation(map(w, :), map(w, :), maxLag);   % nanXcorrelation's ftn is used.
    acf = acf0(maxLag+1:end);
    acfMat1(w, :) = acf;
end

