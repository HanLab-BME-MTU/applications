%[cohorts res] = plotIntensityCohorts(data, varargin) plot average intensities for a range of lifetime cohorts

% Francois Aguet (last modified 08/14/2012)

function [cohorts res] = plotIntensityCohorts(data, varargin)

nCh = numel(data(1).channels);

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addOptional('ch', nCh:-1:1);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('CohortBounds_s', [10 20 40 60 80 100 120]);
ip.addParamValue('ShowVariation', true, @islogical);
ip.addParamValue('Mode', 'percentiles', @(x) any(strcmpi(x, {'SEM', 'percentiles'})));
ip.addParamValue('ShowBackground', false, @islogical);
ip.addParamValue('Rescale', true, @islogical);
ip.addParamValue('RescalingReference', 'med', @(x) any(strcmpi(x, {'max', 'med'})));
ip.addParamValue('ScaleSlaveChannel', true, @islogical);
ip.addParamValue('MaxIntensityThreshold', 0);
ip.addParamValue('LineStyle', '-');
ip.addParamValue('DisplayMode', '');
ip.addParamValue('TrackIndex', []);
ip.addParamValue('Cutoff_f', 5);
ip.addParamValue('Alpha', 0.05);
ip.addParamValue('YTick', []);
ip.addParamValue('RemoveOutliers', false, @islogical);
ip.addParamValue('ShowLegend', false, @islogical);
ip.parse(data, varargin{:});
cohortBounds = ip.Results.CohortBounds_s;


% if no specific channel is selected, all channels are shown
chVec = ip.Results.ch;

mCh = find(strcmp(data(1).source, data(1).channels));
sCh = setdiff(1:nCh, mCh);

nd = numel(data);
kLevel = norminv(1-ip.Results.Alpha/2, 0, 1);

nc = numel(cohortBounds)-1;
b = 5;
framerate = data(1).framerate;

lftData = getLifetimeData(data, 'Overwrite', ip.Results.Overwrite);

% Scale max. intensity distributions
lftFields = {'A', 'sbA', 'ebA', 'sigma_r', 'sbSigma_r', 'ebSigma_r'};
offset = zeros(nCh,nd);
if ip.Results.Rescale && ~isfield(lftData, 'a');
    for c = 1:nCh
        A = arrayfun(@(i) i.A(i.lifetime_s(i.catIdx==1)>=ip.Results.Cutoff_f,:), lftData, 'UniformOutput', false);
        maxA_all = cellfun(@(i) nanmax(i,[],2)', A, 'UniformOutput', false);

        %maxA_all = arrayfun(@(i) nanmax(i.A(:,:,c),[],2), lftData, 'UniformOutput', false);
        [a offset(c,:)] = rescaleEDFs(maxA_all, 'Display', true, 'Reference', ip.Results.RescalingReference, 'FigureName', ['Channel ' num2str(c)]);
        
        % apply scaling
        for i = 1:nd
            for f = 1:numel(lftFields)
                lftData(i).(lftFields{f})(:,:,c) = a(i) * lftData(i).(lftFields{f})(:,:,c);
            end
        end
    end
end
% no need to exclude tracks < cutoff, smallest cohort is [10..19]
if ~isempty(ip.Results.TrackIndex)
    for i = 1:nd
        for f = 1:numel(lftFields)
            lftData(i).(lftFields{f}) = lftData(i).(lftFields{f})(ip.Results.TrackIndex{i},:,:);
        end
        lftData(i).lifetime_s([lftData(i).catIdx]~=1) = [];
        lftData(i).lifetime_s = lftData(i).lifetime_s(ip.Results.TrackIndex{i});
        lftData(i).trackLengths([lftData(i).catIdx]~=1) = [];
        lftData(i).trackLengths = lftData(i).trackLengths(ip.Results.TrackIndex{i});
        lftData(i).catIdx = ones(size(lftData(i).lifetime_s));
    end
end

% test for outliers
if nd>4 && ip.Results.RemoveOutliers
    outlierIdx = [];
    for c = 1:nCh
        maxA_all = arrayfun(@(i) nanmax(i.A(:,:,c),[],2), lftData, 'UniformOutput', false);
        cOut = detectEDFOutliers(maxA_all, offset(c,:), 'FigureName', ['Outliers, channel ' num2str(c)]);
        outlierIdx = [outlierIdx cOut]; %#ok<AGROW>
        if ~isempty(cOut)
            fprintf('Outlier data sets for channel %d:\n', c);
            for i = 1:numel(cOut)
                fprintf('%s\n', getShortPath(data(cOut(i))));
            end
        end
    end
    outlierIdx = unique(outlierIdx);
    if ~isempty(outlierIdx)
        data(outlierIdx) = [];
        lftData(outlierIdx) = [];
        nd = numel(data);
        fprintf('Outliers excluded from intensity cohorts. Indexes: %s\n', num2str(outlierIdx));
    end
end


% loop through data sets, generate cohorts for each
res(1:nd) = struct('interpTracks', [], 'interpSigLevel', []);

% # data points in cohort (including buffer frames)
iLength = arrayfun(@(c) floor(mean(cohortBounds([c c+1]))/framerate) + 2*b, 1:nc);
% time vectors for cohorts
cT = arrayfun(@(i) (-b:i-b-1)*framerate, iLength, 'UniformOutput', false);

for i = 1:nd
    lifetime_s = lftData(i).lifetime_s([lftData(i).catIdx]==1);
    trackLengths = lftData(i).trackLengths([lftData(i).catIdx]==1);
    
    % for intensity threshold in master channel
    maxA = max(lftData(i).A(:,:,mCh), [], 2)';
    
    for ch = 1:nCh % channels
        % interpolate tracks to mean cohort length
        for c = 1:nc % cohorts
            % tracks in current cohort (above threshold)
            cidx = find(cohortBounds(c)<=lifetime_s & lifetime_s<cohortBounds(c+1) & maxA > ip.Results.MaxIntensityThreshold);
            nt = numel(cidx);

            interpTracks = zeros(nt,iLength(c));
            sigma_rMat = zeros(nt,iLength(c));
            cLengths = trackLengths(cidx);
            % loop through track lengths within cohort
            for t = 1:nt
                A = [lftData(i).sbA(cidx(t),:,ch) lftData(i).A(cidx(t),1:cLengths(t),ch) lftData(i).ebA(cidx(t),:,ch)];
                bgr = [lftData(i).sbSigma_r(cidx(t),:,ch) lftData(i).sigma_r(cidx(t),1:cLengths(t),ch) lftData(i).ebSigma_r(cidx(t),:,ch)];
                
                % align to track start
                %w = min(numel(A),iLength);
                %interpTracks(t,1:w) = A(1:w);
                %sigma_r_Ia(t,1:w) = bgr(1:w);
                
                % interpolate to mean length
                xi = linspace(1,cLengths(t)+2*b, iLength(c));
                %interpTracks(t,:) = interp1(1:cLengths(t)+2*b, A, xi, 'cubic');
                interpTracks(t,:) = binterp(A, xi);
                %sigma_rMat(t,:) = interp1(1:cLengths(t)+2*b, bgr, xi, 'cubic');
                sigma_rMat(t,:) = binterp(bgr, xi);
            end

            res(i).interpTracks{ch,c} = interpTracks;
            res(i).interpSigLevel{ch,c} = kLevel*sigma_rMat;
         
            % split as a function of slave channel signal
            if isfield(lftData(i), 'significantSignal')
                sigIdx = lftData(i).significantSignal(ch,lftData(i).catIdx==1)==1;
                res(i).sigIdx{ch,c} = sigIdx(cidx); 
            end
        end
    end
end


fset = loadFigureSettings(ip.Results.DisplayMode);

% Set colormap depending on # channels
cmap = cell(1,nCh);
cv = cell(1,nCh);
if nCh==1
    %cmap = ones(nc,3);
    %cmap(:,1) = (nc:-1:1)/nc;
    %cmap = hsv2rgb(cmap);
    cmap{1} = jet(nc);
    cv{1} = rgb2hsv(cmap{1});
    cv{1}(:,2) = 0.2;
    cv{1} = hsv2rgb(cv{1});
else
    hues = getFluorophoreHues(data(1).markers);
    for ch = 1:nCh
        v = mod(hues(ch)+linspace(-0.1, 0.1, nc)', 1);
        cmap{ch} = hsv2rgb([v ones(nc,1) 0.8*ones(nc,1)]);
        cv{ch} = hsv2rgb([v 0.4*ones(nc,1) ones(nc,1)]);
        %cmap{ch} = repmat(hsv2rgb([hues(ch) 1 0.8]), [nc 1]);
        %cv{ch} = repmat(hsv2rgb([hues(ch) 0.4 1]), [nc 1]);
    end
end

% scale slave channels relative to master (for visualization only)
sf = ones(1,nCh);
if ip.Results.ScaleSlaveChannel% && nd > 1
    for ch = 1:nCh
        iSF = zeros(1,nc);
        for c = 1:nc
            % find largest mean of all cohorts
            M = arrayfun(@(x) mean(x.interpTracks{ch,c},1), res, 'UniformOutput', false);
            M = mean(vertcat(M{:}), 1);
            iSF(c) = max(M);
        end
        sf(ch) = max(iSF);
    end
end



figure(fset.fOpts{:}, 'Name', 'Intensity cohorts');
axes(fset.axOpts{:});
hold on;
A = cell(1,nc);
for ch = chVec
    for c = nc:-1:1
        if nd > 1
            % means for each data set
            AMat = arrayfun(@(x) mean(x.interpTracks{ch,c},1), res, 'UniformOutput', false);
            AMat = vertcat(AMat{:});
            A{ch,c} = nanmean(AMat,1);
            SEM = nanstd(AMat,[],1)/sqrt(nd);
            Amin = A{ch,c} - SEM;
            Aplus = A{ch,c} + SEM;
        else
            % if input is a single data set, show median + percentiles
            M = prctile(res(1).interpTracks{ch,c}, [25 50 75], 1);
            A{ch,c} = M(2,:);
            Amin = M(1,:);
            Aplus = M(3,:);
        end
        if ip.Results.ShowVariation
            fill([cT{c} cT{c}(end:-1:1)], sf(mCh)/sf(ch)*[Amin Aplus(end:-1:1)], cv{ch}(c,:), 'EdgeColor', cmap{ch}(c,:));
        end
        cohorts.t{c} = cT{c};
    end
end
cohorts.A = A;
for ch = chVec
    % Plot mean/median in front
    for c = nc:-1:1
        plot(cT{c}, sf(mCh)/sf(ch)*A{ch,c}, ip.Results.LineStyle, 'Color', cmap{ch}(c,:), 'LineWidth', 1.5);
    end
    
    % Plot signifcance threshold in front
    if ip.Results.ShowBackground && ch==mCh
        % Background level: median of all detections
        if nd>1
            % median background level per cohort for each data set
            medM = arrayfun(@(i) cellfun(@(x) nanmedian(x(:)), i.interpSigLevel(ch,:)) , res, 'UniformOutput', false);
            medM = vertcat(medM{:});
            plot([-10 120], nanmean(medM(:))*[1 1], 'k--', 'LineWidth', 1);
        else
            % median background level per cohort
            medC = cellfun(@(x) median(x(:)), res.interpSigLevel(ch,:));
            plot([-10 120], mean(medC)*[1 1], 'k--', 'LineWidth', 1);
        end
    end
end
set(gca, 'XLim', [-b*framerate-5 cohortBounds(end)], 'XTick', 0:20:200);
if ~isempty(ip.Results.YTick)
    set(gca, 'YTick', ip.Results.YTick, 'YLim', ip.Results.YTick([1 end]));
end
xlabel('Time (s)', fset.lfont{:});
ylabel('Fluo. intensity (A.U.)', fset.lfont{:});


%%
% indiv. figures for cargo+ / cargo-: split based on significance of slave channel
if isfield(res(1), 'sigIdx') && nCh==2
    figure(fset.fOpts{:}, 'Name', 'Intensity cohorts, cargo-positive tracks');
    axes(fset.axOpts{:})
    hold on;
    A = cell(1,nc);
    % plot slave channel first
    ch = 2;
    for c = nc:-1:1
        if nd > 1
            % means for each data set
            AMat = arrayfun(@(x) mean(x.interpTracks{ch,c}(x.sigIdx{2,c},:),1), res, 'UniformOutput', false);
            %AMat = arrayfun(@(x) mean(x.interpTracks{ch,c}(:,:),1), res, 'UniformOutput', false);
            AMat = vertcat(AMat{:});
            A{ch,c} = mean(AMat,1);
            SEM = std(AMat,[],1)/sqrt(nd);
            Amin = A{ch,c} - SEM;
            Aplus = A{ch,c} + SEM;
        else
            % if input is a single data set, show median + percentiles
            M = prctile(res(1).interpTracks{ch,c}(res(1).sigIdx{2,c},:), [25 50 75], 1);
            A{ch,c} = M(2,:);
            Amin = M(1,:);
            Aplus = M(3,:);
        end
        if ip.Results.ShowVariation
            fill([cT{c} cT{c}(end:-1:1)], sf(mCh)/sf(ch)*[Amin Aplus(end:-1:1)], cv{ch}(c,:), 'EdgeColor', cmap{ch}(c,:));
        end
    end
    % plot master channel
    ch = 1;
    for c = nc:-1:1
        if nd > 1
            % means for each data set
            AMat = arrayfun(@(x) mean(x.interpTracks{ch,c}(x.sigIdx{2,c},:),1), res, 'UniformOutput', false);
            AMat = vertcat(AMat{:});
            A{ch,c} = mean(AMat,1);
            SEM = std(AMat,[],1)/sqrt(nd);
            Amin = A{ch,c} - SEM;
            Aplus = A{ch,c} + SEM;
        else
            % if input is a single data set, show median + percentiles
            M = prctile(res(1).interpTracks{ch,c}(res(1).sigIdx{2,c},:), [25 50 75], 1);
            A{ch,c} = M(2,:);
            Amin = M(1,:);
            Aplus = M(3,:);
        end
        if ip.Results.ShowVariation
            fill([cT{c} cT{c}(end:-1:1)], sf(mCh)/sf(ch)*[Amin Aplus(end:-1:1)], cv{ch}(c,:), 'EdgeColor', cmap{ch}(c,:));
        end
    end
    % Plot mean/median in front    
    for ch = [2 1]
        for c = nc:-1:1
            plot(cT{c}, sf(mCh)/sf(ch)*A{ch,c}, ip.Results.LineStyle, 'Color', cmap{ch}(c,:), 'LineWidth', 1.5);
        end
    end
    set(gca, 'XLim', [-b*framerate-5 cohortBounds(end)], 'XTick', 0:20:200);
    if ~isempty(ip.Results.YTick)
        set(gca, 'YTick', ip.Results.YTick, 'YLim', ip.Results.YTick([1 end]));
    end
    xlabel('Time (s)', fset.lfont{:});
    ylabel('Fluo. intensity (A.U.)', fset.lfont{:});
    
    
    %%
    figure(fset.fOpts{:}, 'Name', 'Intensity cohorts, cargo-negative tracks');
    if ip.Results.ShowLegend
        pos = get(gcf, 'Position');
        pos(3) = 10;
        set(gcf, 'Position', pos);
    end
    axes(fset.axOpts{:});
    hold on;
    A = cell(1,nc);
    % plot slave channel first
    ch = 2;
    for c = nc:-1:1
        if nd > 1
            % means for each data set
            AMat = arrayfun(@(x) mean(x.interpTracks{ch,c}(~x.sigIdx{2,c},:),1), res, 'UniformOutput', false);
            AMat = vertcat(AMat{:});
            A{ch,c} = mean(AMat,1);
            SEM = std(AMat,[],1)/sqrt(nd);
            Amin = A{ch,c} - SEM;
            Aplus = A{ch,c} + SEM;
        else
            % if input is a single data set, show median + percentiles
            M = prctile(res(1).interpTracks{ch,c}(~res(1).sigIdx{2,c},:), [25 50 75], 1);
            A{ch,c} = M(2,:);
            Amin = M(1,:);
            Aplus = M(3,:);
        end
        if ip.Results.ShowVariation
            fill([cT{c} cT{c}(end:-1:1)], sf(mCh)/sf(ch)*[Amin Aplus(end:-1:1)], cv{ch}(c,:), 'EdgeColor', cmap{ch}(c,:), 'HandleVisibility', 'off');
        end
    end
    % plot master channel
    ch = 1;
    for c = nc:-1:1
        if nd > 1
            % means for each data set
            AMat = arrayfun(@(x) mean(x.interpTracks{ch,c}(~x.sigIdx{2,c},:),1), res, 'UniformOutput', false);
            AMat = vertcat(AMat{:});
            A{ch,c} = mean(AMat,1);
            SEM = std(AMat,[],1)/sqrt(nd);
            Amin = A{ch,c} - SEM;
            Aplus = A{ch,c} + SEM;
        else
            % if input is a single data set, show median + percentiles
            M = prctile(res(1).interpTracks{ch,c}(~res(1).sigIdx{2,c},:), [25 50 75], 1);
            A{ch,c} = M(2,:);
            Amin = M(1,:);
            Aplus = M(3,:);
        end
        if ip.Results.ShowVariation
            fill([cT{c} cT{c}(end:-1:1)], sf(mCh)/sf(ch)*[Amin Aplus(end:-1:1)], cv{ch}(c,:), 'EdgeColor', cmap{ch}(c,:), 'HandleVisibility', 'off');
        end
    end
    % Plot mean/median in front   
    hp = zeros(1,2*nc);
    for ch = [2 1]
        for c = nc:-1:1
            hp(c + nc*(ch-1)) = plot(cT{c}, sf(mCh)/sf(ch)*A{ch,c}, ip.Results.LineStyle, 'Color', cmap{ch}(c,:), 'LineWidth', 1.5);
        end
    end
    set(gca, 'XLim', [-b*framerate-5 cohortBounds(end)], 'XTick', 0:20:200);
    if ~isempty(ip.Results.YTick)
        set(gca, 'YTick', ip.Results.YTick, 'YLim', ip.Results.YTick([1 end]));
    end
    xlabel('Time (s)', fset.lfont{:});
    ylabel('Fluo. intensity (A.U.)', fset.lfont{:});
    
    if ip.Results.ShowLegend
        cohortLabels = arrayfun(@(i) [' ' num2str(cohortBounds(i)) '-' num2str(cohortBounds(i+1)-framerate) ' s'], 1:nc, 'UniformOutput', false);
        hl = legend(hp, [cohortLabels cohortLabels], 'Location', 'SouthEast');
        set(hl, 'Box', 'off', fset.tfont{:}, 'Position', [8 1.5 1.5 3.5]);
    end
end
