function plotMatchedTracks(data, masterTrack, slaveTracks)


hues = getFluorophoreHues(data.markers);
trackColor = hsv2rgb([hues(2) 1 0.8]);
fillLight = hsv2rgb([hues(2) 0.4 1]);
fillDark = hsv2rgb([hues(2) 0.2 1]);

fillLightBuffer = hsv2rgb([hues(2) 0.4 0.85]);
fillDarkBuffer = hsv2rgb([hues(2) 0.2 0.85]);
kLevel = norminv(1-0.05/2.0, 0, 1); % ~2 std above background


pos = get(0, 'DefaultFigurePosition');
ssize = get(0, 'ScreenSize');
b = min(ssize(4)-pos(4), 400);
pos(2) = pos(2) - b;
pos(4) = pos(4) + b;
figure('Position', pos);

ha(1) = axes('Position', [0.15 0.7 0.8 0.25]);
plotTrack(data, masterTrack, 1, 1, 'Handle', ha(1), 'Legend', 'hide');

ha(2) = axes('Position', [0.15 0.4 0.8 0.25]);
plotTrack(data, masterTrack, 1, 2, 'Handle', ha(2), 'Legend', 'hide');

ha(3) = axes('Position', [0.15 0.1 0.8 0.25]);

for k = 1:length(slaveTracks)
    
    % Plot track
    A = slaveTracks(k).A(1,:);
    c = slaveTracks(k).c(1,:);
    sigma_r = slaveTracks(k).sigma_r(1,:);
    t = (slaveTracks(k).start-1:slaveTracks(k).end-1)*data.framerate;
    
    % alpha = 0.05 level
    fill([t t(end:-1:1)], [c c(end:-1:1)+kLevel*sigma_r(end:-1:1)],...
        fillDark, 'EdgeColor', 'none', 'Parent', ha(3));
    hold(ha(3), 'on');

    % plot amplitude std.
    sigma_a = slaveTracks(k).A_pstd(1,:);
    
    rev = c+A-sigma_a;
    fill([t t(end:-1:1)], [c+A+sigma_a rev(end:-1:1)],...
        fillLight, 'EdgeColor', 'none', 'Parent', ha(3));
    
    % plot track
    ampl = A+c;
    plot(ha(3), t, ampl, '.-', 'Color', trackColor, 'LineWidth', 1);
    
    % plot background level
    plot(ha(3), t, c, '-', 'Color', trackColor);
end
set(ha(3), 'XLim', get(ha(2), 'XLim'), 'YLim', get(ha(2), 'YLim'));

