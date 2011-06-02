% plotFrame(data, tracks, frameIdx, ch, varargin)
%
% INPUTS:   data : data structure
%         tracks : track structure
%       frameIdx : index of the frame
%             ch : channel #
%     {varargin} : optional inputs:
%                      'Visible' : {'on'} | 'off' toggles figure visibility
%                      'Mode' : {'raw'} | 'rgb' | 'mask'
%                      'DynamicRange' : cell array of dynamic ranges for all channels
%                      'Handle' : h, axis handle (for plotting from within GUI)
%                      'VisibleTracks' : {'current'} | 'all'
%                      'ScaleBar' : scale | {[]} display scale bar
%                      'Print' : 'on' | {'off'} generates an EPS in 'data.source/Figures/'

% Francois Aguet, March 16 2011

function plotFrame(data, tracks, frameIdx, ch, varargin)

%======================================
% Parse inputs, set defaults
%======================================
nx = data.imagesize(2);
ny = data.imagesize(1);
psize = data.pixelSize/data.M;
nCh = length(ch);

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addRequired('tracks');
ip.addRequired('frameIdx');
ip.addRequired('ch');
ip.addParamValue('visible', 'on', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.addParamValue('mode', 'raw', @(x) strcmpi(x, 'raw') | strcmpi(x, 'rgb') | strcmpi(x, 'mask'));
ip.addParamValue('print', 'off', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.addParamValue('iRange', cell(1,nCh), @(x) iscell(x));
ip.addParamValue('visibleTracks', 'current', @(x) strcmpi(x, 'current') | strcmpi(x, 'all'));
ip.addParamValue('scaleBar', []);
ip.addParamValue('handle', []);
ip.parse(data, tracks, frameIdx, ch, varargin{:});

if ~isempty(ip.Results.handle)
    ha = ip.Results.handle;
    standalone = false;
else
    h = figure('Visible', 'off', 'PaperPositionMode', 'auto');
    position = get(h, 'Position');
    position(4) = ceil(ny/nx*position(3));
    set(h, 'Position', position, 'Visible', ip.Results.visible);
    ha = axes('Position', [0 0 1 1]);
    standalone = true;
end

xa = 1:nx;
ya = 1:ny;

%======================================
% Plot frame
%======================================
switch ip.Results.mode
    case 'RGB'
        if nCh>3
            error('Max. 3 channels in RGB mode.');
        end
        frame = zeros(ny,nx,3);
        idxRGB = getRGBindex(data.markers);
        for c = 1:nCh
            frame(:,:,idxRGB(c)) = scaleContrast(double(imread(data.framePaths{c}{frameIdx})), ip.Results.iRange{c});
        end
        frame = uint8(frame);
    case 'mask'
        if nCh>1
            error('Mask overlay mode only supports 1 channel.');
        end
        % Display mask only where available
        if (exist(data.maskPaths{frameIdx}, 'file')==2)
            frame = double(imread(data.framePaths{ch}{frameIdx}));
            mask = double(imread(data.maskPaths{frameIdx}));
            frame = rgbOverlay(frame, mask, [1 0 0], ip.Results.iRange{ch});
        else
            frame = double(imread(data.framePaths{ch}{frameIdx}));
        end
    otherwise % grayscale frame
        if nCh>1
            error('Grayscale mode only supports 1 channel.');
        end
        frame = double(imread(data.framePaths{ch}{frameIdx}));
end

%     [sy sx] = size(frame);
%     if sy>sx
%         frame = imrotate(frame, 90);
%     end

imagesc(xa, ya, frame, 'Parent', ha);
colormap(gray(256));
axis(ha, 'image');


%======================================
% Plot tracks
%======================================
if ~isempty(tracks)
    %cmap = jet(data.movieLength);
    cmap = jet(max([tracks.end]-[tracks.start]+1));
    hold(ha, 'on');
    switch ip.Results.visibleTracks
        case 'current'
            idx = find([tracks.start] <= frameIdx & frameIdx <= [tracks.end]);
            for k = idx
                fi = 1:frameIdx-tracks(k).start+1;
                nf = length(tracks(k).t);
                plot(ha, tracks(k).x(ch,1), tracks(k).y(ch,1), '*', 'Color', cmap(nf,:), 'MarkerSize', 5);
                plot(ha, tracks(k).x(ch,fi), tracks(k).y(ch,fi), '-', 'Color', cmap(nf,:));
            end
        case 'all'
            for k = 1:length(tracks)
                nf = length(tracks(k).t);
                plot(ha, tracks(k).x(ch,1), tracks(k).y(ch,1), '*', 'Color', cmap(nf,:), 'MarkerSize', 5);
                plot(ha, tracks(k).x(ch,:), tracks(k).y(ch,:), '-', 'Color', cmap(nf,:));
            end
    end
    hold(ha, 'off');
end

if standalone
    axis(ha, 'off');
end
    
if ~isempty(ip.Results.scaleBar)
    plotScaleBar(ip.Results.scaleBar/psize, 'Label', '5 µm');
end


%======================================
% Print EPS
%======================================
if strcmpi(ip.Results.print, 'on')
    fpath = [data.source 'Figures' filesep];
    if ~(exist(fpath, 'dir')==7)
        mkdir(fpath);
    end
    print(h, '-depsc2', '-loose', [fpath 'frame_' num2str(frameIdx) '_ch' num2str(ch) '.eps']);   
end

if strcmp(ip.Results.visible, 'off')
    close(h);
end
