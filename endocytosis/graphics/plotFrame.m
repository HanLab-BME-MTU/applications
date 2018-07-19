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

% Francois Aguet, March 16 2011 (last modified: 02/06/2012)

function ha = plotFrame(data, tracks, frameIdx, ch, varargin)

%======================================
% Parse inputs, set defaults
%======================================
nx = data.imagesize(2);
ny = data.imagesize(1);
psize = data.pixelSize/data.M;
nCh = length(ch);
mCh = find(strcmp(data.channels, data.source));

if ~isempty(tracks)
    nt = numel(tracks);
else
    nt = 0;
end

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addRequired('tracks');
ip.addRequired('frameIdx');
ip.addRequired('ch', @(x) numel(x)<=3);
ip.addParamValue('Visible', 'on', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.addParamValue('Mode', 'raw', @(x) strcmpi(x, 'raw') | strcmpi(x, 'rgb') | strcmpi(x, 'mask'));
ip.addParamValue('Print', 'off', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.addParamValue('iRange', cell(1,numel(data.channels)), @(x) iscell(x));
ip.addParamValue('DisplayType', 'lifetime', @(x) any(strcmpi(x, {'category', 'lifetime', 'object type', 'random', 'projection'})));
ip.addParamValue('ShowEvents', false, @islogical);
% ip.addParamValue('ShowDetection', false, @islogical); % add optional load fct?
ip.addParamValue('Detection', []);
ip.addParamValue('ShowGaps', true, @islogical);
ip.addParamValue('ScaleBar', []);
ip.addParamValue('ScaleBarLabel', []);
ip.addParamValue('Handle', []);
ip.addParamValue('Colormap', hsv2rgb([rand(nt,1) ones(nt,2)]), @(x) size(x,1)==nt && size(x,2)==3);
ip.parse(data, tracks, frameIdx, ch, varargin{:});

detection = ip.Results.Detection;

if ~isempty(ip.Results.Handle)
    ha = ip.Results.Handle;
    standalone = false;
else
    pos = get(0, 'DefaultFigurePosition');
    pos(4) = ceil(ny/nx*pos(3));
    h = figure('Position', pos, 'Visible', ip.Results.Visible, 'PaperPositionMode', 'auto');
    ha = axes('Position', [0 0 1 1]);
    standalone = true;
    set(h,'DefaultLineLineSmoothing', 'on'); % points are not rendered !!
    set(h,'DefaultPatchLineSmoothing', 'on');
end

frame = zeros(ny,nx,nCh,'uint16');
if ~isfield(data, 'maxProj')
    for c = 1:nCh
        if iscell(data.framePaths{1})
            frame(:,:,c) = imread(data.framePaths{ch(c)}{frameIdx});
        else
            frame(:,:,c) = readtiff(data.framePaths{ch(c)}, frameIdx);
        end
    end
else
    for c = 1:nCh
        frame(:,:,c) = imread(data.maxProj{ch(c)}, 'Index', frameIdx);
    end
end

%======================================
% Plot frame
%======================================
switch ip.Results.Mode
    case 'RGB'
        rframe = zeros(ny,nx,3,'uint8');
        idxRGB = getRGBindex(data.markers);
        for c = 1:nCh
            rframe(:,:,idxRGB(c)) = uint8(scaleContrast(double(frame(:,:,c)), ip.Results.iRange{c}));
        end
        frame = rframe;
        % set channel index to master
        ch = mCh;        
    case 'mask'
        if nCh>1
            error('Mask overlay mode only supports 1 channel.');
        end
        % Display mask only where available
        if ch==mCh && iscell(data.maskPaths) && (exist(data.maskPaths{frameIdx}, 'file')==2)
            mask = double(imread(data.maskPaths{frameIdx}));
            frame = rgbOverlay(double(frame), mask, [1 0 0], ip.Results.iRange{ch});
        elseif ch==mCh && (exist(data.maskPaths, 'file')==2)
            mask = double(readtiff(data.maskPaths, frameIdx));
            frame = rgbOverlay(double(frame), mask, [1 0 0], ip.Results.iRange{ch});
        else
            colormap(gray(256));
        end
    otherwise % grayscale frame
        if nCh>1
            error('Grayscale mode only supports 1 channel.');
        end
        colormap(gray(256));
end

imagesc(frame, 'Parent', ha);
if ~isempty(ip.Results.iRange{ch})
    caxis(ha, ip.Results.iRange{ch});
end
axis(ha, 'image');

if ~isempty(detection) && ~isempty(detection.x) % temporary fix: channel 1 is always selected (assumed master) -> change in call
    hold(ha, 'on');
    isPSF = detection.isPSF(1,:)==1;
    if any(isPSF)
        plot(ha, detection.x(1,isPSF), detection.y(1,isPSF), 'o', 'Color', [0 0.6 0], 'MarkerSize', 8);
    end
    if any(~isPSF)
        plot(ha, detection.x(1,~isPSF), detection.y(1,~isPSF), 'o', 'Color', [0.6 0 0], 'MarkerSize', 8);
    end 
    hold(ha, 'off');
end


%======================================
% Plot tracks
%======================================
if ~isempty(tracks)
    
    trackStarts = [tracks.start];
    trackEnds = [tracks.end];
    
    % visible in current frame
    if strcmpi(ip.Results.DisplayType, 'projection')
        idx = 1:numel(tracks);
    else
        idx = [tracks.start]<=frameIdx & frameIdx<=[tracks.end];
    end
    if sum(idx)>0
        
        np = arrayfun(@(i) numel(i.t), tracks); % points in each track
        nt = numel(tracks);
        
        maxn = max(np);
        X = NaN(maxn, nt);
        Y = NaN(maxn, nt);
        F = NaN(maxn, nt);
        G = zeros(maxn, nt);
        
        for k = 1:nt
            i = 1:np(k);
            X(i,k) = tracks(k).x(ch,:);
            Y(i,k) = tracks(k).y(ch,:);
            F(i,k) = tracks(k).f;
            G(i,k) = tracks(k).gapVect;
        end
        
        if ~strcmpi(ip.Results.DisplayType, 'projection')
            X(F>frameIdx) = NaN;
            Y(F>frameIdx) = NaN;
        end
        
        hold(ha, 'on');
        switch lower(ip.Results.DisplayType)
            case 'lifetime'
                
                % colormap: blue to red, [0..120] seconds
                lifetimes_f = round([tracks.lifetime_s]/data.framerate);
                df = data.movieLength-round(120/data.framerate);
                dcoord = 0.25/df;
                cmap = [jet(round(120/data.framerate)); (0.5:-dcoord:0.25+dcoord)' zeros(df,2)];
                
                M = cmap(lifetimes_f(idx),:);
                set(ha, 'ColorOrder', M);
                plot(ha, X, Y);
               
            case 'category'
                % Categories
                % Ia)  Single tracks
                % Ib)  Single tracks, rejected
                % Ic)  Single tracks cut at beginning or end
                % Id)  Single tracks, persistent
                % IIa) Compound tracks with valid gaps
                % IIb) Compound tracks with invalid gaps
                % IIc) Compound tracks cut at beginning or end
                % IId) Compound tracks, persistent

                idx_Ia = [tracks.catIdx]==1;
                idx_Ib = [tracks.catIdx]==2;
                idx_Ic = [tracks.catIdx]==3;
                idx_Id = [tracks.catIdx]==4;
                idx_IIa = [tracks.catIdx]==5;
                idx_IIb = [tracks.catIdx]==6;
                idx_IIc = [tracks.catIdx]==7;
                idx_IId = [tracks.catIdx]==8;
                
                plot(ha, X(:,idx_Ia), Y(:,idx_Ia), 'Color', [0 1 0]);
                plot(ha, X(:,idx_Ib), Y(:,idx_Ib), 'Color', [1 1 0]);
                plot(ha, X(:,idx_Ic), Y(:,idx_Ic), 'Color', [1 0.5 0]);
                plot(ha, X(:,idx_Id), Y(:,idx_Id), 'Color', [1 0 0]);
                plot(ha, X(:,idx_IIa), Y(:,idx_IIa), 'Color', [0 1 1]);
                plot(ha, X(:,idx_IIb), Y(:,idx_IIb), 'Color', [0 0.5 1]);
                plot(ha, X(:,idx_IIc), Y(:,idx_IIc), 'Color', [0 0 1]);
                plot(ha, X(:,idx_IId), Y(:,idx_IId), 'Color', [0.5 0 1]);
            case 'object type'    
                isCCP = [tracks.isCCP];
                plot(ha, X(:,isCCP), Y(:,isCCP), 'Color', [0 0.8 0]);
                plot(ha, X(:,~isCCP), Y(:,~isCCP), 'Color', [0.8 0 0]);
            case 'random'
                set(ha, 'ColorOrder', ip.Results.Colormap(idx,:));
                plot(ha, X, Y);
            case 'projection'
                
                
                lifetimes_f = round([tracks.lifetime_s]/data.framerate);
                df = data.movieLength-round(120/data.framerate);
                dcoord = 0.25/df;
                cmap = [jet(round(120/data.framerate)); (0.5:-dcoord:0.25+dcoord)' zeros(df,2)];
                
                M = cmap(lifetimes_f(idx),:);
                set(ha, 'ColorOrder', M);
                plot(ha, X, Y);
                
        end
        if ip.Results.ShowEvents
            
            % Births
            bcoord = arrayfun(@(i) [i.x(ch,1) i.y(ch,1)], tracks(trackStarts==frameIdx), 'UniformOutput', false);
            bcoord = vertcat(bcoord{:});
            plot(ha, bcoord(:,1), bcoord(:,2), '*', 'Color', 'g', 'MarkerSize', 8, 'LineWidth', 1);
            
            % Deaths
            dcoord = arrayfun(@(i) [i.x(ch,1) i.y(ch,1)], tracks(trackEnds==frameIdx), 'UniformOutput', false);
            dcoord = vertcat(dcoord{:});
            plot(ha, dcoord(:,1), dcoord(:,2), 'x', 'Color', 'r', 'MarkerSize', 8, 'LineWidth', 1);
        end
        
        if ip.Results.ShowGaps
            x_gap = X(F==frameIdx & G==1);
            y_gap = Y(F==frameIdx & G==1);
            plot(ha, x_gap, y_gap, 'o', 'Color', 'w', 'MarkerSize', 6, 'LineWidth', 1);
        end
    end
    hold(ha, 'off');
end

if standalone
    set(ha, 'TickDir', 'out');
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
    axis(ha, 'off');
    if strcmpi(ip.Results.Mode, 'RGB')==1
        chLabel = '_RGB';
    else
        chLabel = ['_ch' num2str(ch)];
    end
    
    print(h, '-depsc2', '-loose', [fpath 'frame_' num2str(frameIdx) chLabel '_' ip.Results.DisplayType '.eps']);   
end

if strcmp(ip.Results.Visible, 'off')
    close(h);
end
