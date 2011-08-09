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

function ha = plotFrame(data, tracks, frameIdx, ch, varargin)

%======================================
% Parse inputs, set defaults
%======================================
nx = data.imagesize(2);
ny = data.imagesize(1);
nt = length(tracks);
psize = data.pixelSize/data.M;
nCh = length(ch);
mCh = find(strcmp(data.channels, data.source));

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addRequired('tracks');
ip.addRequired('frameIdx');
ip.addRequired('ch');
ip.addParamValue('Visible', 'on', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.addParamValue('Mode', 'raw', @(x) strcmpi(x, 'raw') | strcmpi(x, 'rgb') | strcmpi(x, 'mask'));
ip.addParamValue('Print', 'off', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.addParamValue('iRange', cell(1,nCh), @(x) iscell(x));
ip.addParamValue('DisplayType', 'lifetime', @(x) any(strcmpi(x, {'lifetime', 'category', 'projection', 'all'})));
ip.addParamValue('ShowEvents', 'true', @islogical);
ip.addParamValue('ScaleBar', []);
ip.addParamValue('ScaleBarLabel', []);
ip.addParamValue('Handle', []);
ip.addParamValue('Colormap', hsv2rgb([rand(nt,1) ones(nt,2)]), @(x) size(x,1)==nt && size(x,2)==3);
ip.parse(data, tracks, frameIdx, ch, varargin{:});

if ~isempty(ip.Results.Handle)
    ha = ip.Results.Handle;
    standalone = false;
else
    h = figure('Visible', 'off', 'PaperPositionMode', 'auto');
    position = get(h, 'Position');
    position(4) = ceil(ny/nx*position(3));
    set(h, 'Position', position, 'Visible', ip.Results.Visible);
    ha = axes('Position', [0 0 1 1]);
    standalone = true;
end

%======================================
% Plot frame
%======================================
switch ip.Results.Mode
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
        % set channel index to master
        ch = mCh;        
    case 'mask'
        if nCh>1
            error('Mask overlay mode only supports 1 channel.');
        end
        
        % Display mask only where available
        if ch==mCh && (exist(data.maskPaths{frameIdx}, 'file')==2)
            frame = double(imread(data.framePaths{ch}{frameIdx}));
            mask = double(imread(data.maskPaths{frameIdx}));
            frame = rgbOverlay(frame, mask, [1 0 0], ip.Results.iRange{ch});
        else
            %frame = double(imread(data.framePaths{ch}{frameIdx}));
            frame = uint8(scaleContrast(repmat(double(imread(data.framePaths{ch}{frameIdx})), [1 1 3]), ip.Results.iRange{ch}));
        end
    otherwise % grayscale frame
        if nCh>1
            error('Grayscale mode only supports 1 channel.');
        end
        %frame = double(imread(data.framePaths{ch}{frameIdx}));
        frame = uint8(scaleContrast(repmat(double(imread(data.framePaths{ch}{frameIdx})), [1 1 3]), ip.Results.iRange{ch}));
end

%     [sy sx] = size(frame);
%     if sy>sx
%         frame = imrotate(frame, 90);
%     end

imagesc(frame, 'Parent', ha);
if ~isempty(ip.Results.iRange{ch})
    caxis(ha, ip.Results.iRange{ch});
end
axis(ha, 'image');


%======================================
% Plot tracks
%======================================
if ~isempty(tracks)
    lifetimes_f = [tracks.end]-[tracks.start]+1;
    maxlft = max(lifetimes_f);
    if maxlft>120
        df = maxlft-120;
        dc = 0.25/df;
        cmap = [jet(120); (0.5:-dc:0.25+dc)' zeros(df,2)];
    else
        cmap = jet(maxlft);
    end    
    
    hold(ha, 'on');
    switch lower(ip.Results.DisplayType)
        case 'lifetime'
            % discard any compound tracks
            tracks(arrayfun(@(t) iscell(t.x), tracks)) = [];
            tracks([tracks.valid]~=1) = [];
            
            X = catTrackFields(tracks, data.movieLength, 'x', ch);
            Y = catTrackFields(tracks, data.movieLength, 'y', ch);
            gapMap = catTrackFields(tracks, data.movieLength, 'gapVect', 1)==1;
            
            idx = find([tracks.start] <= frameIdx & frameIdx <= [tracks.end]); 
            if ~isempty(idx)
                M = cmap(lifetimes_f(idx),:);
                set(ha, 'ColorOrder', M);
                plot(ha, X(idx,1:frameIdx)', Y(idx,1:frameIdx)');
                if ip.Results.ShowEvents
                    
                    % Births
                    k = find([tracks(idx).start]==frameIdx);
                    plot(ha, X(idx(k),frameIdx), Y(idx(k),frameIdx), '*', 'Color', 'g', 'MarkerSize', 8, 'LineWidth', 1);
                    
                    % Deaths
                    k = find([tracks(idx).end]==frameIdx);                    
                    plot(ha, X(idx(k),frameIdx), Y(idx(k),frameIdx), 'x', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 1);
                    % Gaps
                    k = find(gapMap(idx,frameIdx))';
                    plot(ha, X(idx(k),frameIdx), Y(idx(k),frameIdx), 'o', 'Color', 'w', 'MarkerSize', 6, 'LineWidth', 1);

                    % Old display: markers have the same color as track
                    % % Births
                    % midx = find([tracks(idx).start]==frameIdx);
                    % for k = midx
                    %     plot(ha, X(idx(k),frameIdx), Y(idx(k),frameIdx), '*', 'Color', M(k,:), 'MarkerSize', 8, 'LineWidth', 1);
                    % end
                    % % Deaths
                    % midx = find([tracks(idx).end]==frameIdx);
                    % for k = midx
                    %     plot(ha, X(idx(k),frameIdx), Y(idx(k),frameIdx), 'x', 'Color', M(k,:), 'MarkerSize', 8, 'LineWidth', 1);
                    % end
                    % % Gaps
                    % midx = find(gapMap(idx,frameIdx))';
                    % for k = midx
                    %     plot(ha, X(idx(k),frameIdx), Y(idx(k),frameIdx), 'o', 'Color', M(k,:), 'MarkerSize', 6, 'LineWidth', 1);
                    % end
                end
            end
        case 'category'
            % plot regular tracks
            cidx = [tracks.valid]==1 & arrayfun(@(t) ~iscell(t.x), tracks);
            X = catTrackFields(tracks(cidx), data.movieLength, 'x', ch);
            Y = catTrackFields(tracks(cidx), data.movieLength, 'y', ch);
            idx = find([tracks(cidx).start] <= frameIdx & frameIdx <= [tracks(cidx).end]);
            plot(ha, [X(idx,1:frameIdx)'; NaN(1, length(idx))],...
                [Y(idx,1:frameIdx)'; NaN(1, length(idx))], 'Color', 'g');
            
            % Persistent tracks
            cidx = arrayfun(@(t) ~iscell(t.x), tracks) & [tracks.status]==3;
            X = catTrackFields(tracks(cidx), data.movieLength, 'x', ch);
            Y = catTrackFields(tracks(cidx), data.movieLength, 'y', ch);
            idx = find([tracks(cidx).start] <= frameIdx & frameIdx <= [tracks(cidx).end]);
            plot(ha, [X(idx,1:frameIdx)'; NaN(1, length(idx))],...
                [Y(idx,1:frameIdx)'; NaN(1, length(idx))], 'Color', 'y');
            
            % Invalid tracks
            cidx = [tracks.valid]==0 & arrayfun(@(t) ~iscell(t.x), tracks) & [tracks.status]~=3;
            X = catTrackFields(tracks(cidx), data.movieLength, 'x', ch);
            Y = catTrackFields(tracks(cidx), data.movieLength, 'y', ch);
            idx = find([tracks(cidx).start] <= frameIdx & frameIdx <= [tracks(cidx).end]);
            plot(ha, [X(idx,1:frameIdx)'; NaN(1, length(idx))],...
                [Y(idx,1:frameIdx)'; NaN(1, length(idx))], 'Color', 'r');
            
            % plot split/merge tracks
            cidx = arrayfun(@(t) iscell(t.x), tracks);
            [X, starts, ends] = catTrackFields(tracks(cidx), data.movieLength, 'x', ch);
            Y = catTrackFields(tracks(cidx), data.movieLength, 'y', ch);
            idx = find(starts <= frameIdx & frameIdx <= ends);
            plot(ha, [X(idx,1:frameIdx)'; NaN(1, length(idx))],...
                [Y(idx,1:frameIdx)'; NaN(1, length(idx))], 'Color', 'c');
            
        case 'projection'
            X = catTrackFields(tracks, data.movieLength, 'x', ch);
            Y = catTrackFields(tracks, data.movieLength, 'y', ch);
            M = cmap(lifetimes_f,:);
            set(ha, 'ColorOrder', M);
            plot(ha, X', Y');
            
        case 'all'
            X = catTrackFields(tracks, data.movieLength, 'x', ch);
            Y = catTrackFields(tracks, data.movieLength, 'y', ch);
            idx = find([tracks.start] <= frameIdx & frameIdx <= [tracks.end]);
            set(ha, 'ColorOrder', ip.Results.Colormap(idx,:));            
            plot(ha, [X(idx,1:frameIdx)'; NaN(1, length(idx))],...
                [Y(idx,1:frameIdx)'; NaN(1, length(idx))]);
            
    end
    hold(ha, 'off');
end

if standalone
    axis(ha, 'off');
end
    
if ~isempty(ip.Results.ScaleBar)
    plotScaleBar(ip.Results.ScaleBar/psize, 'Label', ip.Results.ScaleBarLabel);
end


%======================================
% Print EPS
%======================================
if strcmpi(ip.Results.Print, 'on')
    fpath = [data.source 'Figures' filesep];
    if ~(exist(fpath, 'dir')==7)
        mkdir(fpath);
    end
    print(h, '-depsc2', '-loose', [fpath 'frame_' num2str(frameIdx) '_ch' num2str(ch) '.eps']);   
end

if strcmp(ip.Results.Visible, 'off')
    close(h);
end
