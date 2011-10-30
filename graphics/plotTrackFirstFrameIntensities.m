% Francois Aguet (last modified 10/30/2011)

function plotTrackFirstFrameIntensities(data, tracks)

mCh = strcmp(data.source, data.channels);

idx = [tracks.nSeg]==1 & [tracks.valid]==1;
tracks = tracks(idx);

trackLengths_f = [tracks.end]-[tracks.start]+1;

% groups: 1 2 3 4 5 6-10 11-20 21-40 41-80
b0 = [1 2 3 4 5 6 11 21 41];
b1 = [1 2 3 4 5 10 20 40 80];
ng = numel(b0);

f1mean = zeros(ng,1);
f1std = zeros(ng,1);
f1bgmean = zeros(ng,1);
f1bgstd = zeros(ng,1);
nTracks = zeros(ng,1);
for k = 1:ng
    A0 = arrayfun(@(t) t.x{1}(mCh,1), tracks(trackLengths_f>=b0(k) & trackLengths_f<=b1(k)));
    f1mean(k) = mean(A0);
    f1std(k) = std(A0);
    nTracks(k) = numel(A0);
    
    c0 = arrayfun(@(t) t.sigma_r{1}(mCh,1), tracks(trackLengths_f>=b0(k) & trackLengths_f<=b1(k)));
    f1bgmean(k) = mean(c0);
    f1bgstd(k) = std(c0);
end

kLevel = norminv(1-0.05/2.0, 0, 1); % ~2 std above background

figure('Position', [440 378 760 320], 'PaperPositionMode', 'auto');
barplot2([kLevel*f1bgmean f1mean], [kLevel*f1bgstd f1std], 'ErrorbarPosition', 'both', ...
    'Color', [1 0.76 0.8; 0.7 0.9 1], 'EdgeColor', [0.87 0 0; 0 0.7 1],...
    'GroupDistance', 1, 'BarWidth', 1, 'Angle', 0,...
    'XLabels', {'1', '2', '3', '4', '5', '6-10', '11-20', '21-40', '41-80'},...
    'XLabel', 'Track lifetime (frames)', 'YLabel', '1^{st} frame intensity (A.U.)');

hl = legend('Background', 'Amplitude', 'Location', 'NorthEastOutside');
set(hl, 'Box', 'off');
