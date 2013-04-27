%[cohorts res] = plotIntensityCohorts(data, varargin) plot average intensities for a range of lifetime cohorts

% Francois Aguet (last modified 08/14/2012)

function [cohorts, res] = plotIntensityCohorts(data, varargin)

nCh = numel(data(1).channels);

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addOptional('ch', nCh:-1:1);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('CohortBounds_s', [10 20 40 60 80 100 120]);
ip.addParamValue('ShowVariation', true, @islogical);
ip.addParamValue('FillMode', 'SEM', @(x) any(strcmpi(x, {'SEM', 'pct'})));
ip.addParamValue('FrontLayer', false, @islogical);
ip.addParamValue('ShowBackground', false, @islogical);
ip.addParamValue('Rescale', true, @islogical);
ip.addParamValue('RescalingReference', 'med', @(x) any(strcmpi(x, {'max', 'med'})));
ip.addParamValue('ScaleSlaveChannel', true, @islogical);
ip.addParamValue('ScalingFactor', ones(1,nCh));
ip.addParamValue('MaxIntensityThreshold', 0);
ip.addParamValue('ExcludeVisitors', true, @islogical);
ip.addParamValue('SlaveName', [], @ischar);
ip.addParamValue('ChannelNames', []);
ip.addParamValue('LineStyle', '-');
ip.addParamValue('Hues', []);
ip.addParamValue('DisplayMode', '');
ip.addParamValue('DisplayAll', false);
ip.addParamValue('TrackIndex', []);
ip.addParamValue('Cutoff_f', 5);
ip.addParamValue('Alpha', 0.05);
ip.addParamValue('YTick', []);
ip.addParamValue('YLim', []);
ip.addParamValue('RemoveOutliers', false, @islogical);
ip.addParamValue('ShowLegend', false, @islogical);
ip.addParamValue('ShowPct', true, @islogical);
ip.addParamValue('LftDataName', 'lifetimeData.mat');
ip.parse(data, varargin{:});
cohortBounds = ip.Results.CohortBounds_s;
sf = ip.Results.ScalingFactor;
hues = ip.Results.Hues;

% if no specific channel is selected, all channels are shown
chVec = ip.Results.ch;
mCh = find(strcmp(data(1).source, data(1).channels));

nd = numel(data);
kLevel = norminv(1-ip.Results.Alpha/2.0, 0, 1);

nc = numel(cohortBounds)-1;
b = 5;
framerate = data(1).framerate;

% # data points in cohort (including buffer frames)
iLength = arrayfun(@(c) floor(mean(cohortBounds([c c+1]))/framerate) + 2*b, 1:nc);
% time vectors for cohorts
cT = arrayfun(@(i) (-b:i-b-1)*framerate, iLength, 'UniformOutput', false);

XLim = [-b*framerate-5 cohortBounds(end)];
YLim = ip.Results.YLim;
if isempty(YLim) && ~isempty(ip.Results.YTick)
   YLim = ip.Results.YTick([1 end]); 
end

lftData = getLifetimeData(data, 'Overwrite', ip.Results.Overwrite,...
    'LifetimeData', ip.Results.LftDataName, 'Scale', ip.Results.Rescale, 'Cutoff_f', 5,...
    'ReturnValidOnly', true, 'ExcludeVisitors', ip.Results.ExcludeVisitors);

if ~isempty(ip.Results.TrackIndex)
    lftFields = fieldnames(lftData);
    i = cellfun(@(f) size(lftData(1).(f),1)==numel(lftData(1).lifetime_s), lftFields);
    lftFields = lftFields(i);
    for i = 1:nd
        for f = 1:numel(lftFields)
            lftData(i).(lftFields{f}) = lftData(i).(lftFields{f})(ip.Results.TrackIndex{i},:,:);
        end
    end
end

% loop through data sets, generate cohorts for each
res(1:nd) = struct('interpTracks', [], 'interpSigLevel', []);
cohortBounds(end) = cohortBounds(end)+framerate;
for i = 1:nd
    
    % for intensity threshold in master channel
    maxA = max(lftData(i).A(:,:,mCh), [], 2);
    
    for ch = 1:nCh % channels
        % interpolate tracks to mean cohort length
        for c = 1:nc % cohorts
            % tracks in current cohort (above threshold)
            cidx = find(cohortBounds(c)<=lftData(i).lifetime_s & lftData(i).lifetime_s<cohortBounds(c+1) &...
                maxA > ip.Results.MaxIntensityThreshold);
            nt = numel(cidx);

            interpTracks = zeros(nt,iLength(c));
            sigma_rMat = zeros(nt,iLength(c));
            cLengths = lftData(i).trackLengths(cidx);
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
            if isfield(lftData(i), 'significantMaster')
                sigIdx = lftData(i).significantMaster(:,ch)==1;
                %sigIdx = lftData(i).significantSlave(:,ch)==1;
                %sigIdx = lftData(i).significantMaster(:,ch)==0 & lftData(i).significantSlave(:,ch)==1;
                res(i).sigIdx{ch,c} = sigIdx(cidx); 
            end
        end
    end
end
cohortLabels = arrayfun(@(i) [num2str(cohortBounds(i)) '-' num2str(cohortBounds(i+1)-framerate) 's'], 1:nc, 'Unif', 0);
XTick = (cohortBounds(1:end-1)+[cohortBounds(2:end-1) cohortBounds(end)-framerate])/2;

fset = loadFigureSettings(ip.Results.DisplayMode);

% Set colormap depending on # channels
cmap = cell(1,nCh);
cv = cell(1,nCh);
if nCh==1
    if isempty(hues)
        hues = getFluorophoreHues(data(1).markers);
    end
    %cmap = ones(nc,3);
    %cmap(:,1) = (nc:-1:1)/nc;
    %cmap = hsv2rgb(cmap);
    v = mod(hues(ch)+linspace(-0.1, 0.1, nc)', 1);
    cmap{ch} = hsv2rgb([v ones(nc,1) 0.9*ones(nc,1)]);
    cv{ch} = hsv2rgb([v 0.4*ones(nc,1) ones(nc,1)]);
    
%     cmap{1} = jet(nc);
%     cv{1} = rgb2hsv(cmap{1});
%     cv{1}(:,2) = 0.2;
%     cv{1} = hsv2rgb(cv{1});
else
    if isempty(hues)
        hues = getFluorophoreHues(data(1).markers);
    end
    for ch = 1:nCh
        %v = mod(hues(ch)+linspace(-0.1, 0.1, nc)', 1);
        v = mod(hues(ch)+linspace(-0.05, 0.05, nc)', 1);
        cmap{ch} = hsv2rgb([v ones(nc,1) 0.9*ones(nc,1)]);
        cv{ch} = hsv2rgb([v 0.4*ones(nc,1) ones(nc,1)]);
        %cmap{ch} = repmat(hsv2rgb([hues(ch) 1 0.8]), [nc 1]);
        %cv{ch} = repmat(hsv2rgb([hues(ch) 0.4 1]), [nc 1]);
    end
end

% scale slave channels relative to master (for visualization only)
if ip.Results.ScaleSlaveChannel && nCh>1
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
    sf = sf(mCh)./sf;
end

% output
cohorts.bounds = cohortBounds;


%%
%==================================================
% Plot cohorts
%==================================================

if ~(isfield(res(1), 'sigIdx') && nCh>1)
    
    figure(fset.fOpts{:}, 'Name', 'Intensity cohorts', 'Position', [2 2 8 6]);
    ha(1) = axes(fset.axOpts{:}, 'Position', [1.5 2 6 3.5]);
    hold on;
    A = cell(1,nc);
    hp = zeros(1,nc*nCh);
    for c = nc:-1:1
        for ch = chVec
            if nd > 1
                % means for each data set
                if strcmpi(ip.Results.FillMode, 'SEM')
                    AMat = arrayfun(@(x) nanmean(x.interpTracks{ch,c},1), res, 'UniformOutput', false);
                    AMat = vertcat(AMat{:});
                    A{ch,c} = nanmean(AMat,1);
                    SEM = nanstd(AMat,[],1)/sqrt(nd);
                    Amin = A{ch,c} - SEM;
                    Aplus = A{ch,c} + SEM;
                else
                    %AMat = arrayfun(@(x) prctile(x.interpTracks{ch,c},50,1), res, 'Unif', 0);
                    %AMat = vertcat(AMat{:});
                    %A{ch,c} = mean(AMat,1);
                    %AMat = arrayfun(@(x) prctile(x.interpTracks{ch,c},25,1), res, 'Unif', 0);
                    %AMat = vertcat(AMat{:});
                    %Amin = mean(AMat,1);
                    %AMat = arrayfun(@(x) prctile(x.interpTracks{ch,c},75,1), res, 'Unif', 0);
                    %AMat = vertcat(AMat{:});
                    %Aplus = mean(AMat,1);
                    
                    %AMat = arrayfun(@(x) nanmedian(x.interpTracks{ch,c},1), res, 'UniformOutput', false);
                    %AMat = vertcat(AMat{:});
                    %A{ch,c} = nanmean(AMat,1);
                    
                    %Amin = prctile(AMat,25,1);
                    %Aplus = prctile(AMat,75,1);
                end
                
            else
                % if input is a single data set, show median + percentiles
                M = prctile(res(1).interpTracks{ch,c}, [25 50 75], 1);
                A{ch,c} = M(2,:);
                Amin = M(1,:);
                Aplus = M(3,:);
            end
            if ip.Results.ShowVariation
                hp(c + nc*(ch-1)) = fill([cT{c} cT{c}(end:-1:1)], sf(ch)*[Amin Aplus(end:-1:1)], cv{ch}(c,:), 'EdgeColor', cmap{ch}(c,:));
            end
            if ~ip.Results.FrontLayer
                plot(cT{c}, sf(ch)*A{ch,c}, ip.Results.LineStyle, 'Color', cmap{ch}(c,:), 'LineWidth', 1);
            end
            cohorts.t{c} = cT{c};
            cohorts.Amin{ch,c} = Amin;
            cohorts.Aplus{ch,c} = Aplus;
        end
    end
    cohorts.A = A;
    for ch = chVec
        % Plot mean/median in front
        if ip.Results.FrontLayer
            for c = nc:-1:1
                plot(cT{c}, sf(mCh)/sf(ch)*A{ch,c}, ip.Results.LineStyle, 'Color', cmap{ch}(c,:), 'LineWidth', 1);
            end
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
    if isempty(YLim)
        YLim = get(gca, 'YLim');
    end
    
    if ~isempty(ip.Results.ChannelNames)
        hl = legend(hp(floor(nc/2):nc:end), ip.Results.ChannelNames{:});
        set(hl, 'Box', 'off', 'Position', [6 4.75 1 0.8]);
    end
    
%     if ip.Results.ShowLegend
%         %xlabel(ha(1), 'Time (s)', fset.lfont{:});
%         %xlabel(ha(2), 'Time (s)', fset.lfont{:});
%         %XTick = 0:20:200;
%     else
%         set(gca, 'XTick', XTick, 'XTickLabel', cohortLabels);
%         rotateXTickLabels(gca, 'AdjustFigure', false);
%         xlabel(gca, 'Lifetime cohort', fset.lfont{:});
%     end
    
    
    %xlabel('Time (s)', fset.lfont{:});
    ylabel('Fluo. intensity (A.U.)', fset.lfont{:});
    
%%
% indiv. figures for cargo+ / cargo-: split based on significance of slave channel
else
    if ip.Results.ShowLegend
        fpos = [2 2 17 5.5];
        aposy = 1.5;
    else
        fpos = [2 2 15 6];
        aposy = 2;
    end    
    figure(fset.fOpts{:}, 'Position', fpos, 'Name', 'Intensity cohorts, cargo-positive tracks');
    ha(1) = axes(fset.axOpts{:}, 'Position', [1.5 aposy 6 3.5]);
    hold on;
    A = cell(nCh,nc);
    % plot slave channel first
    
    for c = nc:-1:1
        ch = chVec(2);
        if nd > 1
            % means for each data set
            AMat = arrayfun(@(x) mean(x.interpTracks{ch,c}(x.sigIdx{chVec(2),c},:),1), res, 'UniformOutput', false);
            AMat = vertcat(AMat{:});
            A{ch,c} = nanmean(AMat,1);
            SEM = nanstd(AMat,[],1)/sqrt(nd);
            Amin = A{ch,c} - SEM;
            Aplus = A{ch,c} + SEM;
        else
            % if input is a single data set, show median + percentiles
            M = prctile(res(1).interpTracks{ch,c}(res(1).sigIdx{chVec(2),c},:), [25 50 75], 1);
            A{ch,c} = M(2,:);
            Amin = M(1,:);
            Aplus = M(3,:);
        end
        if ip.Results.ShowVariation
            fill([cT{c} cT{c}(end:-1:1)], sf(ch)*[Amin Aplus(end:-1:1)], cv{ch}(c,:), 'EdgeColor', cmap{ch}(c,:));
        end
        plot(cT{c}, sf(ch)*A{ch,c}, ip.Results.LineStyle, 'Color', cmap{ch}(c,:), 'LineWidth', 1);
        cohorts.t{c} = cT{c};
        cohorts.Amin{ch,c} = Amin;
        cohorts.Aplus{ch,c} = Aplus;
        cohorts.A{ch,c} = A{ch,c};
        %end
        % plot master channel
        ch = chVec(1);
        %for c = nc:-1:1
        if nd > 1
            % means for each data set
            AMat = arrayfun(@(x) mean(x.interpTracks{ch,c}(x.sigIdx{chVec(2),c},:),1), res, 'UniformOutput', false);
            AMat = vertcat(AMat{:});
            A{ch,c} = nanmean(AMat,1);
            SEM = nanstd(AMat,[],1)/sqrt(nd);
            Amin = A{ch,c} - SEM;
            Aplus = A{ch,c} + SEM;
        else
            % if input is a single data set, show median + percentiles
            M = prctile(res(1).interpTracks{ch,c}(res(1).sigIdx{chVec(2),c},:), [25 50 75], 1);
            A{ch,c} = M(2,:);
            Amin = M(1,:);
            Aplus = M(3,:);
        end
        if ip.Results.ShowVariation
            fill([cT{c} cT{c}(end:-1:1)], sf(ch)*[Amin Aplus(end:-1:1)], cv{ch}(c,:), 'EdgeColor', cmap{ch}(c,:));
        end
        plot(cT{c}, sf(ch)*A{ch,c}, ip.Results.LineStyle, 'Color', cmap{ch}(c,:), 'LineWidth', 1);
        cohorts.t{c} = cT{c};
        cohorts.Amin{ch,c} = Amin;
        cohorts.Aplus{ch,c} = Aplus;
        cohorts.A{ch,c} = A{ch,c};
    end
    
    for ch = chVec
        % Plot mean/median in front
        %for c = nc:-1:1
        %    plot(cT{c}, sf(mCh)/sf(ch)*A{ch,c}, ip.Results.LineStyle, 'Color', cmap{ch}(c,:), 'LineWidth', 1);
        %end
        
        % Plot signifcance threshold in front
        if ip.Results.ShowBackground && ch~=mCh
            % Background level: median of all detections
            if nd>1
                % median background level per cohort for each data set
                medM = arrayfun(@(i) cellfun(@(x) nanmedian(x(:)), i.interpSigLevel(ch,:)) , res, 'UniformOutput', false);
                medM = vertcat(medM{:});
                plot([-10 120], sf(ch)*nanmean(medM(:))*[1 1], 'k--', 'LineWidth', 1);
            else
                % median background level per cohort
                medC = cellfun(@(x) median(x(:)), res.interpSigLevel(ch,:));
                plot([-10 120], mean(medC)*[1 1], 'k--', 'LineWidth', 1);
            end
        end
    end
    if isempty(YLim)
        YLim = get(gca, 'YLim');
    end
   
    % pct pos CCPs/cohort
    M = arrayfun(@(r) cellfun(@(i) sum(i)/numel(i)*100, r.sigIdx(chVec(2),:)), res, 'unif', 0);
    % pct of all CCPs
    %M = arrayfun(@(r) cellfun(@(i) sum(i)/numel([r.sigIdx{2,:}])*100, r.sigIdx(2,:)), res, 'unif', 0);
    M = vertcat(M{:});
   
    if ~isempty(ip.Results.SlaveName)
        sv = arrayfun(@(i) 100*sum(i.significantSignal(:,2)==1 & max(i.A(:,:,1),[],2)>ip.Results.MaxIntensityThreshold)...
            /sum(max(i.A(:,:,1),[],2)>ip.Results.MaxIntensityThreshold), lftData);
        text(XLim(1)+0.025*diff(XLim), YLim(2), [ip.Results.SlaveName ' pos. ('  num2str(mean(sv), '%.1f') ' ± ' num2str(std(sv), '%.1f') '% of CCPs)'], fset.sfont{:}, 'VerticalAlignment', 'bottom');
    end
    

   if ip.Results.ShowPct
        pos = get(gcf, 'Position');
        pos(3) = 19;
        set(gcf, 'Position', pos);
    end
    
    ha(2) = axes(fset.axOpts{:}, 'Position', [8.25 aposy 6 3.5]);
    hold on;
    A = cell(1,nc);
    % plot slave channel first
    
    for c = nc:-1:1
        ch = chVec(2);
        if nd > 1
            % means for each data set
            AMat = arrayfun(@(x) mean(x.interpTracks{ch,c}(~x.sigIdx{chVec(2),c},:),1), res, 'UniformOutput', false);
            AMat = vertcat(AMat{:});
            A{ch,c} = nanmean(AMat,1);
            SEM = nanstd(AMat,[],1)/sqrt(nd);
            Amin = A{ch,c} - SEM;
            Aplus = A{ch,c} + SEM;
        else
            % if input is a single data set, show median + percentiles
            M = prctile(res(1).interpTracks{ch,c}(~res(1).sigIdx{chVec(2),c},:), [25 50 75], 1);
            A{ch,c} = M(2,:);
            Amin = M(1,:);
            Aplus = M(3,:);
        end
        if ip.Results.ShowVariation
            fill([cT{c} cT{c}(end:-1:1)], sf(ch)*[Amin Aplus(end:-1:1)], cv{ch}(c,:), 'EdgeColor', cmap{ch}(c,:), 'HandleVisibility', 'off');
        end
        hp(c + nc*(ch-1)) = plot(cT{c}, sf(ch)*A{ch,c}, ip.Results.LineStyle, 'Color', cmap{ch}(c,:), 'LineWidth', 1);
    end
    % plot master channel
    
    for c = nc:-1:1
        ch = chVec(1);
        if nd > 1
            % means for each data set
            AMat = arrayfun(@(x) mean(x.interpTracks{ch,c}(~x.sigIdx{chVec(2),c},:),1), res, 'UniformOutput', false);
            AMat = vertcat(AMat{:});
            A{ch,c} = nanmean(AMat,1);
            SEM = nanstd(AMat,[],1)/sqrt(nd);
            Amin = A{ch,c} - SEM;
            Aplus = A{ch,c} + SEM;
        else
            % if input is a single data set, show median + percentiles
            M = prctile(res(1).interpTracks{ch,c}(~res(1).sigIdx{chVec(2),c},:), [25 50 75], 1);
            A{ch,c} = M(2,:);
            Amin = M(1,:);
            Aplus = M(3,:);
        end
        if ip.Results.ShowVariation
            fill([cT{c} cT{c}(end:-1:1)], sf(ch)*[Amin Aplus(end:-1:1)], cv{ch}(c,:), 'EdgeColor', cmap{ch}(c,:), 'HandleVisibility', 'off');
        end
        hp(c + nc*(ch-1)) = plot(cT{c}, sf(ch)*A{ch,c}, ip.Results.LineStyle, 'Color', cmap{ch}(c,:), 'LineWidth', 1);
    end
    % Plot mean/median in front   
    %hp = zeros(1,2*nc);
    %     for ch = [2 1]
    %         for c = nc:-1:1
    %             hp(c + nc*(ch-1)) = plot(cT{c}, sf(mCh)/sf(ch)*A{ch,c}, ip.Results.LineStyle, 'Color', cmap{ch}(c,:), 'LineWidth', 1);
    %         end
    %     end
    
    %set(gca, 'YTick', [], 'YColor', 'w');
    set(ha(2), 'YTickLabel', []);
    
    
    if ~isempty(ip.Results.SlaveName)
        text(XLim(1)+0.025*diff(XLim), YLim(2), [ip.Results.SlaveName ' neg.'], fset.sfont{:}, 'VerticalAlignment', 'bottom');
    end
    
    if ip.Results.ShowLegend
        cohortLabels = arrayfun(@(i) [' ' num2str(cohortBounds(i)) '-' num2str(cohortBounds(i+1)-framerate) ' s'], 1:nc, 'UniformOutput', false);
        hl = legend(hp, [cohortLabels cohortLabels], 'Location', 'SouthEast');
        set(hl, 'Box', 'off', fset.tfont{:}, 'Position', [6.75+7.65 1.5 1.25 3.5]);
    end
    
    if ip.Results.ShowPct
        axes(fset.axOpts{:}, 'Position', [15.5 2 3 2.5], 'TickLength', fset.TickLength*6/3);
        barplot2(mean(M,1)', std(M,[],1)', 'Angle', 0, 'BarWidth', 1, 'GroupDistance', 1,...
            'FaceColor', 0.8*[1 1 1], 'EdgeColor', 0.4*[1 1 1],...
            'YLim', [0 100], 'LineWidth', 1);
        set(gca, 'FontSize', 8);
        
        h = title(['% ' ip.Results.SlaveName ' pos. CCPs'], fset.sfont{:});
        %h = ylabel('% CCPs/cohort', fset.lfont{:});
        pos = get(h, 'Position');
        %pos(1) = 0.8*pos(1);
        pos(2) = 1.1*pos(2);
        set(h, 'Position', pos);
        set(gca, 'YTick', 0:20:100, 'XTickLabel', cohortLabels);
        rotateXTickLabels(gca, 'AdjustFigure', false);
        xlabel('Lifetime cohort', fset.lfont{:});
    end
    
   
    ylabel(ha(1), 'Fluo. intensity (A.U.)', fset.lfont{:}); 
end

set(ha, 'XLim', XLim, 'XTick', XTick, 'YLim', YLim);
if ~isempty(ip.Results.YTick)
    set(ha, 'YTick', ip.Results.YTick);
end

if ip.Results.ShowLegend
    arrayfun(@(x) xlabel(x, 'Time (s)', fset.lfont{:}), ha);
    %XTick = 0:20:200;
else
    set(ha, 'XTick', XTick, 'XTickLabel', cohortLabels);
    for i = 1:numel(ha)
        rotateXTickLabels(ha(i), 'AdjustFigure', false);
        xlabel(ha(i), 'Lifetime cohort', fset.lfont{:});
    end
end


% if ip.Results.ShowLegend
%     %xlabel(ha(1), 'Time (s)', fset.lfont{:});
%     %xlabel(ha(2), 'Time (s)', fset.lfont{:});
%     %XTick = 0:20:200;
% else
%     set(ha(1), 'XTick', XTick, 'XTickLabel', cohortLabels);
%     rotateXTickLabels(ha(1), 'AdjustFigure', false);
%     xlabel(ha(1), 'Lifetime cohort', fset.lfont{:});
% end


