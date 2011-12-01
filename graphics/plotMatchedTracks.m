function plotMatchedTracks(data, masterTrack, slaveTracks)


hues = getFluorophoreHues(data.markers);
trackColor = hsv2rgb([hues(2) 1 0.8]);
fillLight = hsv2rgb([hues(2) 0.4 1]);
fillDark = hsv2rgb([hues(2) 0.2 1]);

fillLightBuffer = hsv2rgb([hues(2) 0.4 0.85]);
fillDarkBuffer = hsv2rgb([hues(2) 0.2 0.85]);
kLevel = norminv(1-0.05/2.0, 0, 1); % ~2 std above background


% pos = get(0, 'DefaultFigurePosition');
% ssize = get(0, 'ScreenSize');
% b = min(ssize(4)-pos(4), 400);
% pos(2) = pos(2) - b;
% pos(4) = pos(4) + b;
% figure('Position', pos, 'PaperPositionMode', 'auto');

figure('Position', [440 378 560 720], 'PaperPositionMode', 'auto');

ha(1) = axes('Position', [0.15 0.72 0.8 0.27]);
hold on;
plotTrack(data, masterTrack, 1, 'Handle', ha(1), 'DisplayMode', 'print');
YLim{1} = get(ha(1), 'YLim');

ha(2) = axes('Position', [0.15 0.42 0.8 0.27]);
hold on;
plotTrack(data, masterTrack, 2, 'Handle', ha(2), 'DisplayMode', 'print');
YLim{2} = get(ha(2), 'YLim');

% plot track fragments
ha(3) = axes('Position', [0.15 0.08 0.8 0.27]);
hold(ha(3), 'on');
for k = 1:length(slaveTracks)
    
    A = slaveTracks(k).A{1}(1,:);
    c = slaveTracks(k).c{1}(1,:);
    bgcorr = nanmean(c);
    c = c-bgcorr;
    sigma_r = slaveTracks(k).sigma_r{1}(1,:);
    t = (slaveTracks(k).start-1:slaveTracks(k).end-1)*data.framerate;
    
    % alpha = 0.05 level
    fill([t t(end:-1:1)], [c c(end:-1:1)+kLevel*sigma_r(end:-1:1)],...
        fillDark, 'EdgeColor', 'none', 'Parent', ha(3));
    

    gapIdx = arrayfun(@(x,y) x:y, slaveTracks(k).gapStarts{1}, slaveTracks(k).gapEnds{1}, 'UniformOutput', false);
    gapIdx = [gapIdx{:}];
    
    % plot amplitude std.
    sigma_a = slaveTracks(k).A_pstd{1}(1,:);
    
    rev = c+A-sigma_a;
    fill([t t(end:-1:1)], [c+A+sigma_a rev(end:-1:1)],...
        fillLight, 'EdgeColor', 'none', 'Parent', ha(3));
    
    % plot track
    ampl = A+c;
    ampl(gapIdx) = NaN;
    plot(ha(3), t, ampl, '.-', 'Color', trackColor, 'LineWidth', 2, 'MarkerSize', 21);
    
    ampl = A+c;
    ampl(setdiff(gapIdx, 1:length(ampl))) = NaN;
    if ~isempty(gapIdx)
        lh(3) = plot(ha(3), t, ampl, '--', 'Color', trackColor, 'LineWidth', 2);
        lh(4) = plot(ha(3), t(gapIdx), A(gapIdx)+c(gapIdx), 'o', 'Color', trackColor,...
            'MarkerFaceColor', 'w', 'LineWidth', 2, 'MarkerSize', 7);
    end
    
    % plot background level
    plot(ha(3), t, c, '-', 'Color', trackColor, 'LineWidth', 2);
end

set(ha(1:2), 'XTickLabel', []);
set(ha, 'FontName', 'Helvetica', 'FontSize', 20, 'LineWidth', 2, 'TickDir', 'out');
xlabel('Time (s)');
hy = arrayfun(@(x) ylabel(x, 'Intensity (A.U.)'), ha);

set(ha(1), 'YLim', YLim{1});
set(ha(2), 'YLim', YLim{2});

set(ha(3), 'XLim', get(ha(2), 'XLim'), 'YLim', get(ha(2), 'YLim'),...
    'YTick', get(ha(2), 'YTick'), 'Box', 'off');


xpos = zeros(1,3);
for k = 1:3
    pos = get(hy(k), 'Position');
    xpos(k) = pos(1);
end
xpos = min(xpos);
for k = 1:3
    pos = get(hy(k), 'Position');
    pos(1) = xpos;
    set(hy(k), 'Position', pos);
end





