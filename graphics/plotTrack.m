% plotTrack(data, tracks, trackIdx, ch, varargin)
%
% INPUTS:   data : data structure
%          track : track structure
%        trackIdx : index of the track
%             ch : channel #
%     {varargin} : optional inputs:
%                      'Visible' : {'on'} | 'off' toggles figure visibility
%                      'Handle' : h, axis handle (for plotting from within GUI)
%                      'Print' : 'on' | {'off'} generates an EPS in 'data.source/Figures/'

% Francois Aguet, March 9 2011 (split from trackDisplayGUI)

function plotTrack(data, tracks, trackIdx, ch, varargin)

%======================================
% Parse inputs, set defaults
%======================================
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addRequired('tracks', @isstruct);
ip.addRequired('trackIdx', @isnumeric);
ip.addRequired('ch', @isnumeric);
ip.addParamValue('visible', 'on', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.addParamValue('print', 'off', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.addParamValue('handle', []);
ip.addParamValue('Legend', 'show', @(x) strcmpi(x, 'show') | strcmpi(x, 'hide'));
ip.parse(data, tracks, trackIdx, ch, varargin{:});

mCh = find(strcmp(data.channels, data.source));

if ~isempty(ip.Results.handle)
    ha = ip.Results.handle;
    standalone = false;
else
    hfig = figure('Visible', ip.Results.visible);
    ha = axes('Position', [0.15 0.15 0.8 0.8]);
    standalone = true;
end

if length(tracks)>1
    track = tracks(trackIdx);
else
    track = tracks;
end

if isfield(track, 'startBuffer') && ~isempty(track.startBuffer)
    bStart = size(track.startBuffer.A,2);
else
    bStart = 0;
end
if isfield(track, 'endBuffer') && ~isempty(track.endBuffer)
    bEnd = size(track.endBuffer.A,2);
else
    bEnd = 0;
end

hues = getFluorophoreHues(data.markers);
trackColor = hsv2rgb([hues(ch) 1 0.8]);
fillLight = hsv2rgb([hues(ch) 0.4 1]);
fillDark = hsv2rgb([hues(ch) 0.2 1]);

fillLightBuffer = hsv2rgb([hues(ch) 0.4 0.85]);
fillDarkBuffer = hsv2rgb([hues(ch) 0.2 0.85]);

% Significance thresholds
% sigmaT = icdf('normal', 1-alpha/2, 0, 1);
%sigmaL = icdf('normal', 0.95, 0, 1); % weaker, single-tailed
%sigmaH = icdf('normal', 0.99, 0, 1);
kLevel = norminv(1-0.05/2.0, 0, 1); % ~2 std above background


% Plot track
lh = NaN(1,9);

A = track.A(ch,:);
c = track.c(ch,:);
%cStd = track.cStd_mask(ch,:);
%cStd = track.cStd_res(ch,:);
sigma_r = track.sigma_r(ch,:);
t = (track.start-1:track.end-1)*data.framerate;

% alpha = 0.05 level
lh(1) = fill([t t(end:-1:1)], [c c(end:-1:1)+kLevel*sigma_r(end:-1:1)],...
    fillDark, 'EdgeColor', 'none', 'Parent', ha);
hold(ha, 'on');

% alpha = 0.01 level
% fill([t t(end:-1:1)], [c+sigmaL*cStd c(end:-1:1)+sigmaH*cStd(end:-1:1)],...
%     fillDark, 'EdgeColor', 'none', 'Parent', ha);

gapIdx = arrayfun(@(x,y) x:y, track.gapStarts, track.gapEnds, 'UniformOutput', false);
gapIdx = [gapIdx{:}];

% plot amplitude std.
sigma_a = track.A_pstd(ch,:);

rev = c+A-sigma_a;
fill([t t(end:-1:1)], [c+A+sigma_a rev(end:-1:1)],...
    fillLight, 'EdgeColor', 'none', 'Parent', ha);
% plot(ha, t, c+A+sigma_a, '-', 'Color', fillLight);
% plot(ha, t, c+A-sigma_a, '-', 'Color', fillLight);

% plot track
ampl = A+c;
if ch==mCh
    ampl(gapIdx) = NaN;
end
lh(2) = plot(ha, t, ampl, '.-', 'Color', trackColor, 'LineWidth', 1);

% plot gaps separately
if ch==mCh
    ampl = A+c;
    ampl(setdiff(gapIdx, 1:length(ampl))) = NaN;
    if ~isempty(gapIdx)
        lh(3) = plot(ha, t, ampl, '--', 'Color', trackColor, 'LineWidth', 1);
        lh(4) = plot(ha, t(gapIdx), A(gapIdx)+c(gapIdx), 'o', 'Color', trackColor, 'MarkerFaceColor', 'w', 'LineWidth', 1);
    end
end

% plot background level
lh(5) = plot(ha, t, c, '-', 'Color', trackColor);




% Plot left buffer
if isfield(track, 'startBuffer') && ~isempty(track.startBuffer)
    A = [track.startBuffer.A(ch,:) track.A(ch,1)];
    c = [track.startBuffer.c(ch,:) track.c(ch,1)];
    
    sigma_a = [track.startBuffer.A_pstd(ch,:) track.A_pstd(ch,1)];
    sigma_r = [track.startBuffer.sigma_r(ch,:) track.sigma_r(ch,1)];
    t = (track.start-bStart-1:track.start-1)*data.framerate;
     
    fill([t t(end:-1:1)], [c c(end:-1:1)+kLevel*sigma_r(end:-1:1)],...
        fillDarkBuffer, 'EdgeColor', 'none', 'Parent', ha);
    
    rev = c+A-sigma_a;
    fill([t t(end:-1:1)], [c+A+sigma_a rev(end:-1:1)],...
        fillLightBuffer, 'EdgeColor', 'none', 'Parent', ha);
    
    lh(6) = plot(ha, t, A+c, '.--', 'Color', trackColor, 'LineWidth', 1);
    lh(7) = plot(ha, t, c, '--', 'Color', trackColor);
end

% Plot right buffer
if isfield(track, 'endBuffer') && ~isempty(track.endBuffer)
    A = [track.A(ch,end) track.endBuffer.A(ch,:)];
    c = [track.c(ch,end) track.endBuffer.c(ch,:)];
    
    sigma_a = [track.A_pstd(ch,end) track.endBuffer.A_pstd(ch,:)];
    sigma_r = [track.sigma_r(ch,end) track.endBuffer.sigma_r(ch,:)];
    t = (track.end-1:track.end+bEnd-1)*data.framerate;
    
    fill([t t(end:-1:1)], [c c(end:-1:1)+kLevel*sigma_r(end:-1:1)],...
        fillDarkBuffer, 'EdgeColor', 'none', 'Parent', ha);

    rev = c+A-sigma_a;
    fill([t t(end:-1:1)], [c+A+sigma_a rev(end:-1:1)],...
        fillLightBuffer, 'EdgeColor', 'none', 'Parent', ha);
    
    lh(8) = plot(ha, t, A+c, '.--', 'Color', trackColor, 'LineWidth', 1);
    lh(9) = plot(ha, t, c, '--', 'Color', trackColor);
end


l = legend(lh([2 5 1]), ['Amplitude ch. ' num2str(ch)], ['Background ch. ' num2str(ch)], '\alpha = 0.95 level', 'Location', 'NorthEast');
legend(l, ip.Results.Legend);


tlength = track.end+bEnd - track.start-bStart + 1;
set(ha, 'XLim', ([track.start-bStart-0.1*tlength track.end+bEnd+0.1*tlength]-1)*data.framerate);
box off;



if standalone
    tfont = {'FontName', 'Helvetica', 'FontSize', 14};
    sfont = {'FontName', 'Helvetica', 'FontSize', 18};
    lfont = {'FontName', 'Helvetica', 'FontSize', 22};
    
    set(l, tfont{:});
    
    set(gca, 'LineWidth', 1.5, sfont{:});
    xlabel('Time (s)', lfont{:})
    ylabel('Intensity (A.U.)', lfont{:});

    for k = lh([2 3 5 6:9])
        if ~isnan(k)
            set(k, 'LineWidth', 2);
        end
    end
    
    for k = lh([2 6 8])
        if ~isnan(k)
            set(k, 'MarkerSize', 21);
        end
    end
    
    if ~isnan(lh(4))
        set(lh(4), 'MarkerSize', 7, 'LineWidth', 2);
    end
end

if strcmpi(ip.Results.print, 'on')
    fpath = [data.source 'Figures' filesep];
    if ~(exist(fpath, 'dir')==7)
        mkdir(fpath);
    end
    print(hfig, '-depsc2', '-r300', [fpath 'track_' num2str(trackIdx) '_ch' num2str(ch) '.eps']);
end

if strcmp(ip.Results.visible, 'off')
    close(hfig);
end

