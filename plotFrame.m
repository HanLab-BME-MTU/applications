% plotFrame(data, tracks, frameIdx, ch, varargin)
%
% INPUTS:   data : data structure
%          track : track structure
%        trackID : index of the track
%             ch : channel #
%     {varargin} : optional inputs:
%                      'Visible' : {'on'} | 'off' toggles figure visibility
%                      'Units' : sets the axis units. Default: pixels
%                      'Mode' : {'raw'} | 'rgb' | 'mask'
%                      'DynamicRange' : cell array of dynamic ranges for all channels
%                      'Handle' : h, axis handle (for plotting from within GUI)
%                      'VisibleTracks' : {'current'} | 'all'
%                      'ScaleBar' : scale | {[]} display scale bar
%                      'Print' : 'on' | {'off'} generates an EPS in 'data.source/Figures/'

% Francois Aguet, March 16 2011

function plotFrame(data, tracks, frameIdx, ch, varargin)

nCh = length(ch);

%======================================
% Check inputs
%======================================
if mod(length(varargin),2)~=0
    error('Optional arguments need to be entered as pairs.');
end

%======================================
% Parse inputs, set defaults
%======================================
nx = data.imagesize(2);
ny = data.imagesize(1);
psize = data.pixelSize/data.M;

idx = find(strcmpi(varargin, 'Visible'));
if ~isempty(idx)
    visible = varargin{idx+1};
else
    visible = 'on';
end

idx = find(strcmpi(varargin, 'Mode'));
if ~isempty(idx)
    mode = varargin{idx+1};
else
    mode = 'raw';
end

idx = find(strcmpi(varargin, 'DynamicRange'));
if ~isempty(idx)
    iRange = varargin{idx+1};
else
    iRange = cell(1,nCh);
end

idx = find(strcmpi(varargin, 'Handle'));
if ~isempty(idx)
    ha = varargin{idx+1};
    standalone = false;
else
    h = figure('Visible', 'off', 'PaperPositionMode', 'auto');
    position = get(h, 'Position');
    position(4) = ceil(ny/nx*position(3));
    set(h, 'Position', position, 'Visible', visible);
    ha = axes('Position', [0 0 1 1]);
    standalone = true;
end

idx = find(strcmpi(varargin, 'Units'));
if ~isempty(idx)
    unit = varargin{idx+1};
    xa = (0:nx-1)*psize*unit;
    ya = (0:ny-1)*psize*unit;
else
    xa = 1:nx;
    ya = 1:ny;
end

idx = find(strcmpi(varargin, 'VisibleTracks'));
if ~isempty(idx)
    visibleTracks = varargin{idx+1};
else
    visibleTracks = 'current';
end

idx = find(strcmpi(varargin, 'Print'));
if ~isempty(idx) && strcmpi(varargin{idx+1}, 'on')
    printEPS = true;
else
    printEPS = false;
end

idx = find(strcmpi(varargin, 'ScaleBar'));
if ~isempty(idx)
    scaleBar = varargin{idx+1};
else
    scaleBar = [];
end


%======================================
% Plot frame
%======================================
switch mode
    case 'RGB'
        if nCh>3
            error('Max. 3 channels in RGB mode.');
        end
        frame = zeros(ny,nx,3);
        idxRGB = assignColorsHSV(data.markers);
        for c = 1:nCh
            frame(:,:,idxRGB(c)) = scaleContrast(double(imread(data.framePaths{c}{frameIdx})), iRange{c});
        end
        frame = uint8(frame);
    case 'mask'
        if nCh>1
            error('Mask overlay mode only supports 1 channel.');
        end
        % Display mask only where available
        if (exist([data.channels{ch} 'Detection' filesep 'detectionResults.mat'], 'file')==2)
            overlayColor = [1 0 0];
            frame = scaleContrast(double(imread(data.framePaths{ch}{frameIdx})), iRange{ch});
            [chR chG chB] = deal(frame);
            mask = double(imread(data.maskPaths{frameIdx}));
            maskIdx = mask~=0;
            chR(maskIdx) = chR(maskIdx)*overlayColor(1);
            chG(maskIdx) = chG(maskIdx)*overlayColor(2);
            chB(maskIdx) = chB(maskIdx)*overlayColor(3);
            frame = uint8(cat(3, chR, chG, chB));
        else
            frame = double(imread(data.framePaths{ch}{frameIdx}));
        end
    otherwise % grayscale frame
        if nCh>1
            error('Grayscale mode only supports 1 channel.');
        end
        frame = double(imread(data.framePaths{ch}{frameIdx}));
end

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
    switch visibleTracks
        case 'current'
            idx = find([tracks.start] <= frameIdx & frameIdx <= [tracks.end]);
            for k = idx
                fi = 1:frameIdx-tracks(k).start+1;
                nf = length(tracks(k).x);
                plot(ha, tracks(k).x(1), tracks(k).y(1), '*', 'Color', cmap(nf,:), 'MarkerSize', 5);
                plot(ha, tracks(k).x(fi), tracks(k).y(fi), '-', 'Color', cmap(nf,:));
            end
        case 'all'
            for k = 1:length(tracks)
                nf = length(tracks(k).x);
                plot(ha, tracks(k).x(1), tracks(k).y(1), '*', 'Color', cmap(nf,:), 'MarkerSize', 5);
                plot(ha, tracks(k).x, tracks(k).y, '-', 'Color', cmap(nf,:));
            end
    end
    hold(ha, 'off');
end

if standalone
    axis(ha, 'off');
end

if ~isempty(scaleBar)
    plotScaleBar(scaleBar/psize, 'Label', '5 µm');
end


%======================================
% Print EPS
%======================================
if printEPS
    fpath = [data.source 'Figures' filesep];
    if ~(exist(fpath, 'dir')==7)
        mkdir(fpath);
    end
    print(h, '-depsc2', '-loose', [fpath 'frame_' num2str(frameIdx) '_ch' num2str(ch) '.eps']);   
end

if strcmp(visible, 'off')
    close(h);
end


function idxRGB = assignColorsHSV(markers)

hue = cellfun(@(x) rgb2hsv(wavelength2rgb(name2wavelength(x))), markers, 'UniformOutput', false);
hue = vertcat(hue{:});
hue = hue(:,1); % retain only 'h' from hsv

[hue, sortIdx] = sort(hue, 'descend');

hueRef = [0 1/3 2/3]; % RGB
N = length(markers);
switch N
    case 1
        J = (hue-hueRef).^2;
        idxRGB = find(J==min(J));
    case {2,3}
        idxRGB = sortIdx;
    otherwise
        error('Max. 3 channels for RGB display.');
end
