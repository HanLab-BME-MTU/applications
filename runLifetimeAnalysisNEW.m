function [res] = runLifetimeAnalysisNEW(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(x) isstruct(x) && numel(unique([data.framerate]))==1);
ip.addParamValue('Display', 'on', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.addParamValue('FileName', 'trackAnalysis.mat', @ischar);
ip.addParamValue('Type', 'all', @ischar);
ip.addParamValue('Cutoff', 4, @isscalar);
ip.addParamValue('Print', false, @islogical);
ip.addParamValue('Buffer', 5);
ip.parse(data, varargin{:});
nd = length(data);
framerate = data(1).framerate;
buffer = ip.Results.Buffer;

% Extend all to max. movie length, in case of mismatch
Nmax = max([data.movieLength])-2;

cutoff_f = ip.Results.Cutoff; %REDUNDANT
cutoff_s = ip.Results.Cutoff * framerate;

% generate lifetime histograms
fprintf('Lifetime analysis:   0%%');
for k = 1:nd
    
    % load tracks
    load([data(k).source 'Tracking' filesep ip.Results.FileName]);
    
    % Apply cut-off
    idx = [tracks.lifetime_s]>=cutoff_s;
    tracks = tracks(idx);
    
    % load/create cell mask, discard tracks that fall into background
    mpath = [data(k).source 'Detection' filesep 'cellmask.tif'];
    if (exist(mpath, 'file')==2)
        mask = logical(imread(mpath));
    else
        mask = logical(getCellMask(data(k)));
    end
    x = round(arrayfun(@(tr) nanmean(tr.x(1,:)), tracks));
    y = round(arrayfun(@(tr) nanmean(tr.y(1,:)), tracks));
    
    % exclude tracks in background
    [ny,nx] = size(mask);
    idx = sub2ind([ny nx], y, x);
    tracks = tracks(mask(idx)==1);
    
    
    lifetimes_s = [tracks.lifetime_s];
    
    %====================
    % Track statistics
    %====================
    % Categories
    % Ia)  Single tracks with valid gaps
    % Ib)  Single tracks with invalid gaps
    % Ic)  Single tracks cut at beginning or end
    % Id)  Single tracks, persistent
    % IIa) Compound tracks with valid gaps
    % IIb) Compound tracks with invalid gaps
    % IIc) Compound tracks cut at beginning or end
    % IId) Compound tracks, persistent
    validGaps = arrayfun(@(t) max([t.gapStatus 4]), tracks)==4;
    singleIdx = [tracks.nSeg]==1;
    vis = [tracks.visibility];
    
    idx_Ia = singleIdx & validGaps & vis==1;
    idx_Ib = singleIdx & ~validGaps & vis==1;
    idx_IIa = ~singleIdx & validGaps & vis==1;
    
    v = [sum(idx_Ia);
        sum(idx_Ib);
        sum(singleIdx & vis==2);
        sum(singleIdx & vis==3);
        sum(idx_IIa);
        sum(~singleIdx & ~validGaps & vis==1);
        sum(~singleIdx & vis==2);
        sum(~singleIdx & vis==3)];
    
    if sum(v) ~= numel(tracks)
        error('Track classification error');
    end
    res.trackClassStats{k} = v/numel(tracks);
    
    %====================
    % Histogram etc.
    %====================
    % longest observable lifetime (in frames): N = movieLength-2*buffer
    N = data(k).movieLength-2*buffer;
    t = (cutoff_f:N)*framerate;
    lftHist_Iab = hist(lifetimes_s(idx_Ia | idx_Ib), t);
    lftHist_Ia = hist(lifetimes_s(idx_Ia), t);
    lftHist_Ib = hist(lifetimes_s(idx_Ib), t);
    lftHist_IIa = hist(lifetimes_s(idx_IIa), t);
    
    % apply correction
    % relative probabilities:
    % P(obs. lifetime==1) = N
    % P(obs. lifetime==N) = 1
    % => weighting:
    w = N./(N-cutoff_f+1:-1:1);
    lftHist_Iab = lftHist_Iab .* w;
    lftHist_Ia = lftHist_Ia .* w;
    lftHist_Ib = lftHist_Ib .* w;
    lftHist_IIa = lftHist_IIa .* w;
    
    % Pad with trailing zeros
    if N<Nmax
        lftHist_Iab = [lftHist_Iab zeros(1,Nmax-N)]; %#ok<AGROW>
        lftHist_Ia = [lftHist_Ia zeros(1,Nmax-N)]; %#ok<AGROW>
        lftHist_Ib = [lftHist_Ib zeros(1,Nmax-N)]; %#ok<AGROW>
        lftHist_IIa = [lftHist_IIa zeros(1,Nmax-N)]; %#ok<AGROW>
    end
    
    % Normalization
    lftHist_Iab = lftHist_Iab / sum(lftHist_Iab) / framerate;
    lftHist_Ia = lftHist_Ia / sum(lftHist_Ia) / framerate;
    lftHist_Ib = lftHist_Ib / sum(lftHist_Ib) / framerate;
    lftHist_IIa = lftHist_IIa / sum(lftHist_IIa) / framerate;
    
    res.lftHist_Iab{k} = lftHist_Iab;
    res.lftHist_Ia{k} = lftHist_Ia;
    res.lftHist_Ib{k} = lftHist_Ib;
    res.lftHist_IIa{k} = lftHist_IIa;
    
    samples = lifetimes_s(idx_Ia);
    res.samples{k} = samples;
    res.nSamples(k) = numel(samples);
    
    % birth/death statistics
    startsPerFrame_all = hist([tracks.start], 1:data(k).movieLength);
    %startsPerFrame_Ia = hist([tracks(idx_Ia).start], 1:data(k).movieLength);
    
    %====================
    % Initiation density
    %====================
    
    % Cell area
    px = data(k).pixelSize / data(k).M; % pixels size in object space
    res.area(k) = sum(mask(:)) * px^2 / 1e-12; % in µm^2
    spf = startsPerFrame_all(4:end-cutoff_f);
    
    madFactor = 1/norminv(0.75, 0, 1);
    res.init_um_min(k,:) = [mean(spf); madFactor*mad(spf, 1)]/data(k).framerate*60/res.area(k);
    
    %====================
    % Gap statistics
    %====================
    
    binEdges = [0:20:120 data(k).movieLength-data(k).framerate];
    nb = length(binEdges)-1;
    gapsPerTrack_Ia = zeros(1,nb);
    gapsPerTrack_Ib = zeros(1,nb);
    gapsPerTrack_IIa = zeros(1,nb);
    for b = 1:nb
        tidx = binEdges(b)<=lifetimes_s & lifetimes_s<binEdges(b+1);
        gapsPerTrack_Ia(b) = mean(arrayfun(@(i) sum(i.gapVect), tracks(idx_Ia & tidx)));
        gapsPerTrack_Ib(b) = mean(arrayfun(@(i) sum(i.gapVect), tracks(idx_Ib & tidx)));
        gapsPerTrack_IIa(b) = mean(arrayfun(@(i) sum(i.gapVect), tracks(idx_IIa & tidx)));
    end
    res.gapsPerTrack_Ia{k} = gapsPerTrack_Ia;
    res.gapsPerTrack_Ib{k} = gapsPerTrack_Ib;
    res.gapsPerTrack_IIa{k} = gapsPerTrack_IIa;
    fprintf('\b\b\b\b%3d%%', round(100*k/(nd)));
end
fprintf('\n');

%-------------------------
% Mean histogram
%-------------------------
t_hist = (cutoff_f:Nmax)*framerate;

% Class percentages
v = mean([res.trackClassStats{:}],2);
v_std = std([res.trackClassStats{:}],[],2);

% meanHist_Iab =  mean(vertcat(res.lftHist_Iab{:}),1);
meanHist_Ia =  mean(vertcat(res.lftHist_Ia{:}),1);
meanHist_Ib =  mean(vertcat(res.lftHist_Ib{:}),1);
meanHist_IIa = mean(vertcat(res.lftHist_IIa{:}),1);



%-------------------------
% Assign output
%-------------------------
res.t = t_hist;
% res.meanHist_Ia = meanHist_Ia;
res.source = {data.source};

if strcmpi(ip.Results.Display, 'on')
    fset = loadFigureSettings();
    
    % classes
    hf(1) = plotTrackClasses(v, v_std, 'FaceColor', fset.cfTrackClasses, 'EdgeColor', fset.ceTrackClasses);
    
    % mean histograms (main classes)
    hf(2) = figure;
    hold on;
    hp(3) = plot(t_hist, meanHist_IIa, '.-', 'Color', 0.6*[1 1 1], 'LineWidth', 2, 'MarkerSize', 16);
    hp(2) = plot(t_hist, meanHist_Ib, '.-', 'Color', hsv2rgb([0 1 0.8]), 'LineWidth', 2, 'MarkerSize', 16);
    hp(1) = plot(t_hist, meanHist_Ia, '.-', 'Color', 'k', 'LineWidth', 2, 'MarkerSize', 16);
    axis([0 min(120, t_hist(end)) 0 0.05]);
    set(gca, 'LineWidth', 2, fset.sfont{:}, 'Layer', 'top');
    xlabel('Lifetime (s)', fset.lfont{:});
    ylabel('Frequency', fset.lfont{:});
    hl = legend(hp, 'Single tracks', 'Single tracks, rej. gaps', 'Comp. track', 'Location', 'NorthEast');
    set(hl, 'Box', 'off', fset.tfont{:});
    
    % gap statistics
    ce = fset.ceTrackClasses([1 2 5],:);
    cf = fset.cfTrackClasses([1 2 5],:);
    xlabels = arrayfun(@(b) [num2str(binEdges(b)) '-' num2str(binEdges(b+1)) ' s'], 1:numel(binEdges)-1, 'UniformOutput', false);
    
    
    M_Ia = vertcat(res.gapsPerTrack_Ia{:});
    M_Ib = vertcat(res.gapsPerTrack_Ib{:});
    M_IIa = vertcat(res.gapsPerTrack_IIa{:});
    
    M = [mean(M_Ia,1); mean(M_Ib,1); mean(M_IIa,1)]';
    S = [std(M_Ia,[],1); std(M_Ib,[],1); std(M_IIa,[],1)]';
    
    hf(3) = figure; barplot2(M, S, 'FaceColor', cf, 'EdgeColor', ce,...
        'XLabels', xlabels, 'XLabel', 'Lifetime cohort', 'YLabel', 'gaps/track');
    
end

if ip.Results.Print
    
    fpath = cell(1,nd);
    for k = 1:nd
        [~,fpath{k}] = getCellDir(data(k));
    end
    fpath = unique(fpath);
    if numel(fpath)>1
        fprintf('Figures could not be printed.');
    else
        fpath = [fpath{1} 'Figures' filesep];
        [~,~] = mkdir(fpath);
        print(hf(1), '-depsc2', [fpath 'trackClassDistribution.eps']);
        print(hf(2), '-depsc2', [fpath 'meanLftHist_classes.eps']);
        print(hf(3), '-depsc2', [fpath 'gapStatistics.eps']);
    end
end
