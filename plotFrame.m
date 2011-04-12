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
ip.addParamValue('print', false, @islogical);
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
    set(h, 'Position', position, 'Visible', visible);
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
    
if ~isempty(ip.Results.scaleBar)
    plotScaleBar(ip.Results.scaleBar/psize, 'Label', '5 µm');
end


%======================================
% Print EPS
%======================================
if ip.Results.print
    fpath = [data.source 'Figures' filesep];
    if ~(exist(fpath, 'dir')==7)
        mkdir(fpath);
    end
    print(h, '-depsc2', '-loose', [fpath 'frame_' num2str(frameIdx) '_ch' num2str(ch) '.eps']);   
end

if strcmp(ip.Results.visible, 'off')
    close(h);
end
