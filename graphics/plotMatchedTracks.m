function ha = plotMatchedTracks(data, masterTrack, slaveTracks, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('DisplayMode', 'print', @(x) any(strcmpi(x, {'print', 'screen'})));
ip.addParamValue('MarkerSizes', [6 1.5 0.5]);
ip.addParamValue('YTick',[]);
ip.addParamValue('YLim', []);
ip.addParamValue('XTick',[]);
ip.addParamValue('XLim',[]);
ip.parse(varargin{:});
markerSizes = ip.Results.MarkerSizes;
lineWidth = 0.5;
YLim = ip.Results.YLim;
XLim = ip.Results.XLim;
YTick = ip.Results.YTick;
XTick = ip.Results.XTick;

dt = data.framerate;
if isempty(XTick)
    XTick = 0:20:data.movieLength*dt+20;
end

fset = loadFigureSettings(ip.Results.DisplayMode);
hues = getFluorophoreHues(data.markers);

ah = 1.5;
aw = 3;
if isempty(XLim)
    XLim = [-7*dt masterTrack.t(end)-masterTrack.t(1)+7*dt];
else
    wref = fset.axPos(3); % reference
    width = masterTrack.t(end)-masterTrack.t(1) +7*dt + 7*dt;
    xscale = width/diff(XLim);
    aw = aw*xscale;
    if aw>2
        tickLength = fset.TickLength*wref/aw;
    else
        tickLength = fset.TickLength*wref/2;
    end
    fset.axOpts = [fset.axOpts, 'TickLength', tickLength];
    XLim = [-7*dt masterTrack.t(end)-masterTrack.t(1)+7*dt];
end



figure(fset.fOpts{:}, 'Position', [6 6 8 10]);

% Plot 'master'
ha(1) = axes(fset.axOpts{:}, 'Position', [1.5 5 aw ah]);
hold on;
plotTrack(data, masterTrack, 1, 'Handle', ha(1), 'DisplayMode', 'print', 'MarkerSizes', markerSizes, 'LineWidth', lineWidth);



% Plot 'slave' (2nd channel)
ha(2) = axes(fset.axOpts{:}, 'Position', [1.5 3.25 aw ah]);
hold on;
plotTrack(data, masterTrack, 2, 'Handle', ha(2), 'DisplayMode', 'print', 'MarkerSizes', markerSizes, 'LineWidth', lineWidth);

if isempty(YLim)
    YLim{1} = get(ha(1), 'YLim');
    YLim{2} = get(ha(2), 'YLim');
end

% plot slaveTracks fragments
ha(3) = axes(fset.axOpts{:}, 'Position', [1.5 1.5 aw ah]);
hold(ha(3), 'on');
for k = 1:length(slaveTracks)
    plotTrack(data, slaveTracks(k), 1, 'Handle', ha(3), 'DisplayMode', 'print',...
        'Time', masterTrack.t(1), 'PlotBuffers', false, 'Hues', hues(2), 'MarkerSizes', markerSizes, 'LineWidth', lineWidth);
end

set(ha(1:2), 'XTickLabel', []);
set(ha, 'XLim', XLim, 'XTick', XTick);
% xlabel('Time (s)');

set(ha(1), 'YLim', YLim{1}, 'YTick', YTick{1});
set(ha(2), 'YLim', YLim{2}, 'YTick', YTick{2});
set(ha(3), 'XLim', XLim, 'YLim', YLim{2},...
    'YTick', YTick{2}, 'Box', 'off');

% hy(2) = ylabel(ha(2), 'Intensity (A.U.)', fset.lfont{:});
% hy([1 3]) = arrayfun(@(x) ylabel(x, ''), ha([1 3]));
% hy = arrayfun(@(x) ylabel(x, 'Intensity (A.U.)'), ha);

% xpos = zeros(1,3);
% for k = 1:3
%     pos = get(hy(k), 'Position');
%     xpos(k) = pos(1);
% end
% xpos = min(xpos);
% for k = 1:3
%     pos = get(hy(k), 'Position');
%     pos(1) = xpos;
%     set(hy(k), 'Position', pos);
% end
