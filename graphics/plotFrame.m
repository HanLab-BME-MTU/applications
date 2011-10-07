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

function ha = plotFrame(data, trackInfo, frameIdx, ch, varargin)

%======================================
% Parse inputs, set defaults
%======================================
nx = data.imagesize(2);
ny = data.imagesize(1);
psize = data.pixelSize/data.M;
nCh = length(ch);
mCh = find(strcmp(data.channels, data.source));

% if trackMatrix is a track structure, convert to matrix;
% this makes selection/visualization of a subset of tracks possible
if isfield(trackInfo, 'A')
    tracks = trackInfo;
    clear('trackInfo');
    trackInfo.x = catTrackFields(tracks, data.movieLength, 'x', mCh);
    trackInfo.y = catTrackFields(tracks, data.movieLength, 'y', mCh);
    [gapMap, segStarts, segEnds, seg2trackIndex, track2segIndex] = catTrackFields(tracks, data.movieLength, 'gapVect', mCh);
    gapMap = gapMap==1;
    trackInfo.gapMap = gapMap;
    trackInfo.segStarts = segStarts;
    trackInfo.segEnds = segEnds;
    trackInfo.seg2trackIndex = seg2trackIndex;
    trackInfo.track2segIndex = track2segIndex;
    trackInfo.nSeg = [tracks.nSeg];
    trackInfo.status = [tracks.status];
    trackInfo.valid = [tracks.valid];
    trackInfo.lifetimes_f = [tracks.end]-[tracks.start]+1;     
end

if ~isempty(trackInfo)
    trackIdx = unique(trackInfo.seg2trackIndex);
    nt = numel(trackIdx);
else
    nt = 0;
end

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addRequired('trackMatrix');
ip.addRequired('frameIdx');
ip.addRequired('ch');
ip.addParamValue('Visible', 'on', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.addParamValue('Mode', 'raw', @(x) strcmpi(x, 'raw') | strcmpi(x, 'rgb') | strcmpi(x, 'mask'));
ip.addParamValue('Print', 'off', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.addParamValue('iRange', cell(1,nCh), @(x) iscell(x));
ip.addParamValue('DisplayType', 'lifetime', @(x) any(strcmpi(x, {'lifetime', 'category', 'projection', 'all'})));
ip.addParamValue('ShowEvents', false, @islogical);
ip.addParamValue('ShowDetection', false, @islogical);
ip.addParamValue('ShowGaps', true, @islogical);
ip.addParamValue('ScaleBar', []);
ip.addParamValue('ScaleBarLabel', []);
ip.addParamValue('Handle', []);
ip.addParamValue('Colormap', hsv2rgb([rand(nt,1) ones(nt,2)]), @(x) size(x,1)==nt && size(x,2)==3);
ip.parse(data, trackInfo, frameIdx, ch, varargin{:});

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

if ip.Results.ShowDetection && ch==mCh
    load([data.source 'Detection' filesep 'detection_v2.mat']);
    hold(ha, 'on');
    plot(ha, frameInfo(frameIdx).x(ch,:), frameInfo(frameIdx).y(ch,:), 'ro', 'MarkerSize', 8);
    hold(ha, 'off');
end


%======================================
% Plot tracks
%======================================
if ~isempty(trackInfo)
     
    
    hold(ha, 'on');
    fi = [1:frameIdx frameIdx];

    switch lower(ip.Results.DisplayType)
        case 'lifetime'
            
            maxlft = max(trackInfo.lifetimes_f);
            if maxlft>120
                df = maxlft-120;
                dc = 0.25/df;
                cmap = [jet(120); (0.5:-dc:0.25+dc)' zeros(df,2)];
            else
                cmap = jet(maxlft);
            end
            
            % single-segment tracks visible in this frame
            idx = find([trackInfo.segStarts] <= frameIdx & frameIdx <= [trackInfo.segEnds] &...
                trackInfo.nSeg(trackInfo.seg2trackIndex)==1);
            if ~isempty(idx)
                M = cmap(trackInfo.lifetimes_f(trackInfo.seg2trackIndex(idx)),:);
                set(ha, 'ColorOrder', M);
                plot(ha, trackInfo.x(idx,fi)', trackInfo.y(idx,fi)');
                if ip.Results.ShowEvents
                    
                    % Births
                    plot(ha, trackInfo.x(trackInfo.segStarts==frameIdx,frameIdx),...
                        trackInfo.y(trackInfo.segStarts==frameIdx,frameIdx), '*', 'Color', 'g', 'MarkerSize', 8, 'LineWidth', 1);
                    
                    % Deaths
                    plot(ha, trackInfo.x(trackInfo.segEnds==frameIdx,frameIdx),...
                        trackInfo.y(trackInfo.segEnds==frameIdx,frameIdx), 'x', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 1);

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
                
                if ip.Results.ShowGaps
                    k = find(trackInfo.gapMap(idx,frameIdx))';
                    plot(ha, trackInfo.x(idx(k),frameIdx), trackInfo.y(idx(k),frameIdx), 'o', 'Color', 'w', 'MarkerSize', 6, 'LineWidth', 1);
                end
            end
        case 'category'
            
            
            % plot regular tracks (single-segment, valid)
            idx = [trackInfo.segStarts] <= frameIdx & frameIdx <= [trackInfo.segEnds] &...
                trackInfo.nSeg(trackInfo.seg2trackIndex)==1 &...
                trackInfo.status(trackInfo.seg2trackIndex)==1 &...
                trackInfo.valid(trackInfo.seg2trackIndex)==1;
            if sum(idx)~=0
                plot(ha, trackInfo.x(idx,fi)', trackInfo.y(idx,fi)', 'Color', 'g');
            end
            
            % Persistent tracks
            idx = trackInfo.status(trackInfo.seg2trackIndex)==3;
            if sum(idx)~=0
                plot(ha, trackInfo.x(idx,fi)',trackInfo.y(idx,fi)', 'Color', 'y');
            end
            
            % Invalid tracks (single-segment)
            idx = [trackInfo.segStarts] <= frameIdx & frameIdx <= [trackInfo.segEnds] &...
                trackInfo.nSeg(trackInfo.seg2trackIndex)==1 &...
                trackInfo.status(trackInfo.seg2trackIndex)~=3 &...
                trackInfo.valid(trackInfo.seg2trackIndex)==0;
            if sum(idx)~=0
                plot(ha, trackInfo.x(idx,fi)', trackInfo.y(idx,fi)', 'Color', 'r');
            end
            
            % plot split/merge tracks
            idx = [trackInfo.segStarts] <= frameIdx & frameIdx <= [trackInfo.segEnds] &...
                trackInfo.nSeg(trackInfo.seg2trackIndex)>1;
            if sum(idx)~=0
                plot(ha, trackInfo.x(idx,fi)', trackInfo.y(idx,fi)', 'Color', 'c');
            end
        case 'all'
            idx = [trackInfo.segStarts] <= frameIdx & frameIdx <= [trackInfo.segEnds];
            
            set(ha, 'ColorOrder', ip.Results.Colormap(trackInfo.seg2trackIndex(idx),:));
            plot(ha, trackInfo.x(idx,fi)', trackInfo.y(idx,fi)');
            
            if ip.Results.ShowEvents
                % Births
                k = trackInfo.segStarts==frameIdx;
                plot(ha, trackInfo.x(k,frameIdx), trackInfo.y(k,frameIdx), '*', 'Color', 'g', 'MarkerSize', 8, 'LineWidth', 1);
                
                % Deaths
                k = trackInfo.segEnds==frameIdx;
                plot(ha, trackInfo.x(k,frameIdx), trackInfo.y(k,frameIdx), 'x', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 1);
            end
            
            if ip.Results.ShowGaps
                idx = idx & gapMap(:,frameIdx);
                plot(ha, trackInfo.x(idx,frameIdx), trackInfo.y(idx,frameIdx), 'o', 'Color', 'w', 'MarkerSize', 6, 'LineWidth', 1);
            end
        case 'projection'
            maxlft = max(trackInfo.lifetimes_f);
            if maxlft>120
                df = maxlft-120;
                dc = 0.25/df;
                cmap = [jet(120); (0.5:-dc:0.25+dc)' zeros(df,2)];
            else
                cmap = jet(maxlft);
            end
            
            set(ha, 'ColorOrder', cmap(trackInfo.lifetimes_f,:));
            plot(ha, trackInfo.x', trackInfo.y');
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
