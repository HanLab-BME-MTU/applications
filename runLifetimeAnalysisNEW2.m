function [res] = runLifetimeAnalysisNEW2(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(x) isstruct(x) && numel(unique([data.framerate]))==1);
ip.addParamValue('Display', 'on', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.addParamValue('FileName', 'trackAnalysis2.mat', @ischar);
ip.addParamValue('Type', 'all', @ischar);
ip.addParamValue('Cutoff', 4, @isscalar);
% ip.addParamValue('Tracks', []);
ip.addParamValue('Print', false, @islogical);
ip.parse(data, varargin{:});
nd = length(data);

framerate = data(1).framerate;

% for k = 1:nd
%     if isempty(ip.Results.Tracks)
%         data(k).tracks = loadTracks(data(k), 'Cutoff', ip.Results.Cutoff, 'Type', ip.Results.Type,...
%             'FileName', ip.Results.FileName);
%     else
%         data(k).tracks = ip.Results.Tracks;
%     end
%     data(k).lifetimes_s = ;
% end

% Extend all to max. movie length, in case of mismatch
Nmax = max([data.movieLength])-2;

cutoff = ip.Results.Cutoff * framerate;

% generate lifetime histograms
for k = 1:nd

    % load tracks
    load([data.source 'Tracking' filesep ip.Results.FileName]);
    
    % Apply cut-off
    idx = [tracks.lifetime_s]>=cutoff;
    tracks = tracks(idx);    
    
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
    res.trackClassDist{k} = v/numel(tracks);
    
    %if sum(data(k).trackStats)~=length(tracks)
    %    error('Selection error');
    %end
        
    %====================
    % Histogram etc.
    %====================
    N = data(k).movieLength-2;
    t = (1:N)*framerate;
    lftHist = hist(lifetimes_s(idx_Ia), t);
    lftHist_Ib = hist(lifetimes_s(idx_Ib), t);
    lftHist_ms = hist(lifetimes_s(idx_IIa), t);
    
    %lftHist = hist(data(k).lifetimes_s([tracks.nSeg]==1 & [tracks.status]==1), t);
    %lftHist_ms = hist(data(k).lifetimes_s([tracks.nSeg]>1 & [tracks.status]==1), t);
    
    % apply correction
    % longest observable lifetime: N = movieLength-2
    % relative probabilities:
    % P(obs. lifetime==1) = N
    % P(obs. lifetime==N) = 1
    % => weighting:
    lftHist = lftHist .* N./(N:-1:1);
    lftHist_Ib = lftHist_Ib .* N./(N:-1:1);
    lftHist_ms = lftHist_ms .* N./(N:-1:1);
    
    % Pad with trailing zeros
    if N<Nmax
        lftHist = [lftHist zeros(1,Nmax-N)]; %#ok<AGROW>
        lftHist_Ib = [lftHist_Ib zeros(1,Nmax-N)]; %#ok<AGROW>
        lftHist_ms = [lftHist_ms zeros(1,Nmax-N)]; %#ok<AGROW>
    end
    
    % Exclude short tracks from histogram
    cutoff_f = ip.Results.Cutoff; %REDUNDANT
    lftHist = lftHist(cutoff_f:end);
    lftHist_Ib = lftHist_Ib(cutoff_f:end);
    lftHist_ms = lftHist_ms(cutoff_f:end);
    % Normalization
    lftHist = lftHist / sum(framerate*lftHist);
    lftHist_Ib = lftHist_Ib / sum(framerate*lftHist_Ib);
    lftHist_ms = lftHist_ms / sum(framerate*lftHist_ms);
    
    res.lftHist{k} = lftHist;
    res.lftHist_Ib{k} = lftHist_Ib;
    res.lftHist_ms{k} = lftHist_ms;
    
    samples = lifetimes_s(idx_Ia);
    res.samples{k} = samples;
    res.nSamples(k) = numel(samples);
    
%     % birth/death statistics
%     starts = [tracks.start];
%     data(k).startsPerFrame = arrayfun(@(f) sum(starts==f), 1:data(k).movieLength);
%     ends = [tracks.end];
%     data(k).endsPerFrame = arrayfun(@(f) sum(ends==f), 1:data(k).movieLength);
% 
%     % idx = [tracks.valid]==1;
%     % starts = [tracks(idx).start];
%     % starts = arrayfun(@(f) sum(starts==f), 1:data(k).movieLength);
%     % ends = [tracks(idx).end];
%     % ends = arrayfun(@(f) sum(ends==f), 1:data(k).movieLength);
%     % data(k).startsPerFrame_valid = starts;
%     % data(k).endsPerFrame_valid = ends;
%     %
%     % idx = [tracks.type]==1 & [tracks.valid]==0 & [tracks.status]==1;
%     % starts = [tracks(idx).start];
%     % starts = arrayfun(@(f) sum(starts==f), 1:data(k).movieLength);
%     % ends = [tracks(idx).end];
%     % ends = arrayfun(@(f) sum(ends==f), 1:data(k).movieLength);
%     % data(k).startsPerFrame_invalid = starts;
%     % data(k).endsPerFrame_invalid = ends;
%     %
%     % idx = [tracks.type]==2 & [tracks.status]==1;
%     % starts = [tracks(idx).start];
%     % starts = arrayfun(@(f) sum(starts==f), 1:data(k).movieLength);
%     % ends = [tracks(idx).end];
%     % ends = arrayfun(@(f) sum(ends==f), 1:data(k).movieLength);
%     % data(k).startsPerFrame_ms = starts;
%     % data(k).endsPerFrame_ms = ends;
%     
%     
%     
%     %====================
%     % Initiation density
%     %====================
%     
%     % Estimate cell area
%     frame = double(imread(data(k).framePaths{1}{1}));
%     bf = bilateralFilter(scaleContrast(frame,[],[0 255]), 5, 100);
%     %figure; imagesc(frame); axis image; colormap(gray(256));
%     %figure; imagesc(bf); axis image; colormap(gray(256));
% 
%     try
%         T = thresholdFluorescenceImage(bf);
%         mask = bf>T;
%     catch
%         mask = ones(size(bf));
%     end
%         
%     %figure; imagesc(mask); axis image; colormap(gray(256));
%     px = data(k).pixelSize / data(k).M;
%     data(k).area = sum(mask(:)) * px^2 / 1e-12; % in µm^2
%     % CC = bwconncomp(mask,8);
%     % csize = cellfun(@numel, CC.PixelIdxList);
%     % CC.PixelIdxList = CC.PixelIdxList(csize == max(csize));
%     % CC.NumObjects = 1;
%     % mask = labelmatrix(CC);
%     %figure; imagesc(mask); axis image; colormap(gray(256)); colorbar;
%     
%     %spf = data(k).startsPerFrame(2:end-cutoff_f+1);
%     %figure; plot(data(k).startsPerFrame)
%     
%     spf = data(k).startsPerFrame(10:end-9);
%     meanInit_f = mean(spf);
%     madFactor = 1/norminv(0.75, 0, 1);
%     %stdInit_f = std(spf);
%     stdInit_f = madFactor*mad(spf, 1);
%     data(k).init_mum_min = [meanInit_f; stdInit_f]/data(k).framerate*60/data(k).area;    
% 
%     %====================
%     % Gap statistics
%     %====================
%     idxValid = [tracks.valid]==1; 
%     idxInvalid = [tracks.nSeg]==1 & [tracks.valid]==0 & [tracks.status]==1;
%     idxMS = [tracks.nSeg]>1 & [tracks.status]==1;
%     
%     lifetimes_s = [tracks.lifetime_s];
%     %lifetimes_f = lifetimes_s/data(k).framerate;
%     bins = [0:20:100 data(k).movieLength];
%     nb = length(bins)-1;
%     
%     %lifetimes_s = [tracks(idxValid).lifetime_s];
%     binIdx = arrayfun(@(b) idxValid & bins(b)<=lifetimes_s & lifetimes_s<bins(b+1), 1:nb, 'UniformOutput', false);
%     for b = 1:nb
%         data(k).gapsPerTrack_valid(b) = mean(arrayfun(@(t) sum([t.gapVect{:}]), tracks(binIdx{b})));
%     end
%     
%     %lifetimes_s = [tracks(idxInvalid).lifetime_s];
%     binIdx = arrayfun(@(b) idxInvalid & bins(b)<=lifetimes_s & lifetimes_s<bins(b+1), 1:nb, 'UniformOutput', false);
%     for b = 1:nb
%         data(k).gapsPerTrack_invalid(b) = mean(arrayfun(@(t) sum([t.gapVect{:}]), tracks(binIdx{b})));
%     end
%     
%     %lifetimes_s = [tracks.lifetime_s];
%     binIdx = arrayfun(@(b) idxMS & bins(b)<=lifetimes_s & lifetimes_s<bins(b+1), 1:nb, 'UniformOutput', false);
%     %binIdx{1}
%     for b = 1:nb
%         data(k).gapsPerTrack_MS(b) = mean(arrayfun(@(t) sum([t.gapVect{:}]), tracks(binIdx{b})));
%     end
%     data(k).bins = bins;
    
end


%-------------------------
% Mean histogram
%-------------------------
t_hist = (cutoff_f:Nmax)*framerate;

% Class percentages
v = mean([res.trackClassDist{:}],2);
v_std = std([res.trackClassDist{:}],[],2);
% plotTrackClasses(v, v_std);



M = vertcat(res.lftHist{:});
meanHist = mean(M,1);
% histSEM = std(M,[],1) / sqrt(length(data));

meanHist_ms = mean(vertcat(res.lftHist_ms{:}),1);
meanHist_Ib = mean(vertcat(res.lftHist_Ib{:}),1);


%-------------------------
% Assign output
%-------------------------
res.t = t_hist;
res.meanHist = meanHist;
res.source = {data.source};


if strcmpi(ip.Results.Display, 'on')
    fset = loadFigureSettings();

    figure;
    hold on;
    %figure('Position', [440 358 560 420]);
    %figure('Position', [440 358 3*dLeft+sum(axesWidths) 370], 'PaperPositionMode', 'auto',...
    %    'InvertHardCopy', 'off', 'Color', [1 1 1], 'Units', 'pixels');
    %axes('Units', 'pixels', 'Position', [dLeft 72 axesWidths(1) 240]);
       
    %     fill([t_hist t_hist(end:-1:1)], [meanHist-histSEM meanHist(end:-1:1)+histSEM(end:-1:1)],...
    %         [1 1 1]*0.7, 'EdgeColor', 'none');
    
    hp(3) = plot(t_hist, meanHist_ms, '.-', 'Color', 0.6*[1 1 1], 'LineWidth', 2, 'MarkerSize', 16);
    hp(2) = plot(t_hist, meanHist_Ib, '.-', 'Color', [0.8 0 0], 'LineWidth', 2, 'MarkerSize', 16);
    hp(1) = plot(t_hist, meanHist, '.-', 'Color', 'k', 'LineWidth', 2, 'MarkerSize', 16);

    axis([0 min(120, t_hist(end)) 0 0.05]);
    set(gca, 'LineWidth', 2, fset.sfont{:}, 'Layer', 'top');
    xlabel('Lifetime (s)', fset.lfont{:});
    ylabel('Frequency', fset.lfont{:});
    
    hl = legend(hp, 'Single tracks', 'Single tracks, rej. gaps', 'Comp. track', 'Location', 'NorthEast');
    set(hl, 'Box', 'off', fset.tfont{:});
    %pos = get(hl, 'Position')
    %pos(2) = 1.5*pos(2);
    %set(hl, 'Position', pos);

end
    
    %     dLeft = 70;
    %
    %     axesWidths = (nbars+1)*25;
    %
    %
    %     insetWidths = (cellfun(@(x) length(x), aVect)+1)*15;
    %
    %
    %
    %     h = bar(x{1}, [steps{1}' model{1}']);
    %
    %     set(h(1), 'BarWidth', 1.5, 'FaceColor', cf1, 'EdgeColor', ce1, 'LineWidth', 1.5);
    %     set(h(2), 'BarWidth', 1.5, 'FaceColor', cf2, 'EdgeColor', ce2, 'LineWidth', 1.5);
    %     axis([-1 x{1}(end)+1 0 ymax])
    %     set(gca, 'LineWidth', 1.5, sfont{:}, 'YTick', 0:0.1:ymax, 'TickDir', 'out', 'Layer', 'top', 'TickLength', [0.01 0]);
    %     xlabel('# EGFP', lfont{:})
    %     ylabel('Frequency', lfont{:});
    %     box off
    %     % text(-0.5, 0.47, ['\gamma = ' num2str(sub_est, '%.2f')], sfont{:});
    %
    %     YLim = get(gca, 'YLim');
    %     h = title('Step 1',lfont{:});
    %     pos = get(h, 'Position');
    %     pos(2) = YLim(2)+(YLim(2)/100)*7;
    %     set(h, 'Position', pos);
    %
    %
    %     % inset
    %     axes('Units', 'pixels', 'Position', [dLeft+axesWidths(1)-insetWidths(1)+offset 240 insetWidths(1) 64]);
    %     h = bar(1:length(aVect{1}), aVect{1});
    %     set(h, 'BarWidth', 0.6, 'FaceColor', cf3, 'EdgeColor', ce3, 'LineWidth', 1.5);
    %     box off;
    %     ya = 0:0.2:1.0;
    %     yal = ['0' arrayfun(@(k) num2str(k, '%.1f'), ya(2:end), 'UniformOutput', false)];
    %     set(gca, 'YTick', ya, 'YTickLabel', yal, 'XTickLabel', nVect{1}/unit, 'LineWidth', 1.5, tfont{:}, 'Layer', 'top', 'TickLength', [0 0]);
    %     axis([0.3 length(aVect{1})+0.7 0 1.05]);
    %     xlabel(['# ' label]);
    %     ylabel('Frequency');
    %
    %
    
    
    
    
%     % Inset
%     hi = axes('Units', 'pixels', 'Position', [360 300 140 80]);
%     
%     M = mean(vertcat(data.trackStats),1);
%     S = std(vertcat(data.trackStats),[],1);
%     
%     barplot2(M, S,...
%         'XLabels', {'valid', 'invalid', 'merge/split', 'persistent', 'cut'}, 'YLabel', '% tracks',...
%         'BarWidth', 0.8, 'FaceColor', 0.8*[1 1 1], 'EdgeColor', 0.6*[1 1 1],...
%         'Handle', hi, 'AdjustFigure', false, 'LabelFontSize', 16, 'AxisFontSize', 14);
%     set(hi, 'YTick', 0:0.2:1, 'YTickLabel', ['0' arrayfun(@(t) num2str(t, '%.1f'), 0.2:0.2:1, 'UniformOutput', false)]);
%     
    
%     hf(2) = figure;    
%     M = [mean(vertcat(data.gapsPerTrack_valid),1);...
%          mean(vertcat(data.gapsPerTrack_invalid),1);...
%          mean(vertcat(data.gapsPerTrack_MS),1)];
%     S = [std(vertcat(data.gapsPerTrack_valid),[],1);...
%          std(vertcat(data.gapsPerTrack_invalid),[],1);...
%          std(vertcat(data.gapsPerTrack_MS),[],1)];
%      
%     xlabels = arrayfun(@(b) [num2str(data(1).bins(b)) '-' num2str(data(1).bins(b+1)-data(1).framerate) ' s'], 1:length(data(1).bins)-1, 'UniformOutput', false);
%     barplot2(M', S', 'FaceColor', [0 0.8 0; 0.8 0 0; 0.6 0.6 0.6], 'XLabels', xlabels, 'YLabel', '# gaps/track');
%     legend('Valid', 'Invalid', 'Merge/split', 'Location', 'NorthWest');
    
    
%     if ip.Results.Print
%         fpath = unique(cellfun(@(c) getParentDir(c),...
%             arrayfun(@(d) getCellDir(d), data, 'UniformOutput', false),...
%             'UniformOutput', false));
%         if numel(fpath)>1
%             fprintf('Figures could not be printed.');
%         else
%             fpath = [fpath{1} 'Figures' filesep];
%             [~,~] = mkdir(fpath);
%             print(hf(1), '-depsc2', [fpath 'lifetimeHist.eps']);
%             print(hf(2), '-depsc2', [fpath 'gapStatistics.eps']);
%         end
%     end
    
    %figure;
    %hold on;
    %plot(data.startsPerFrame(2:end-cutoff_f+1), 'g');
    %plot(data.endsPerFrame, 'r');
    %plot(data.startsPerFrame_valid, 'g');
    %plot(data.endsPerFrame_valid, 'r');
    %plot(data.startsPerFrame_invalid, 'b');
    %plot(data.endsPerFrame_invalid, 'c');
    %plot(data.startsPerFrame_ms, 'm');
    %plot(data.endsPerFrame_ms, 'k');
    

