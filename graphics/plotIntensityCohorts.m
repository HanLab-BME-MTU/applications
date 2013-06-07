%[cohorts res] = plotIntensityCohorts(data, varargin) plot average intensities for a range of lifetime cohorts

% Francois Aguet (last modified 04/30/2013)

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
ip.addParamValue('SlaveName', [], @iscell);
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
ip.addParamValue('AvgFun', @nanmean, @(x) isa(x, 'function_handle'));
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
                res(i).sigIdx{c}(:,ch) = sigIdx(cidx);
            else
                res(i).sigIdx{c}(:,ch) = ones(numel(cidx),1);
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
    if nCh==2
        hb = 0.1;
    else
        hb = 0.05;
    end
    for ch = 1:nCh
        v = mod(hues(ch)+linspace(-hb, hb, nc)', 1);
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



%==================================================
% Plot cohorts
%==================================================
switch nCh
    case 1
        ah = 1;
        na = 1;
        sigCombIdx = [];
    case 2
        na = 2;
        ah = 1;
        sigCombIdx = [1 0]';
        
        pct = zeros(nd,2);
        for i = 1:nd
            s =  lftData(i).significantMaster;
            pct(i,:) = sum([s(:,2) ~s(:,2)],1)/size(s,1);
        end
        meanPct = mean(pct,1);
        stdPct = std(pct,[],1);
    case 3
        na = 4;
        ah = 2;
        sigCombIdx = [1 1; 1 0; 0 1; 0 0];
        
        pct = zeros(nd,4);
        for i = 1:nd
            s = lftData(i).significantMaster;
            vidx = max(lftData(i).A(:,:,mCh),[],2) > ip.Results.MaxIntensityThreshold;
            s = s(vidx,:);
            pct(i,:) = sum([s(:,2)&s(:,3) s(:,2)&~s(:,3) ~s(:,2)&s(:,3) ~s(:,2)&~s(:,3)],1)/size(s,1);
        end
        meanPct = mean(pct,1);
        stdPct = std(pct,[],1);
end

if ~isempty(sigCombIdx)
    SlaveName = ip.Results.SlaveName;
    if isempty(SlaveName)
        SlaveName = data(1).markers(2:nCh);
    end
    tmp = sigCombIdx;
    tmp(tmp==1) = '+';
    tmp(tmp==0) = '-';
    atext = cell(1,na);
    switch nCh
        case 2
            for a = 1:na
                atext{a} = [tmp(a,1) SlaveName{1} ': ' num2str(meanPct(a)*100, '%.1f') '±' num2str(stdPct(a)*100, '%.1f') '%'];
            end
        case 3
            for a = 1:na
                atext{a} = [SlaveName{1} tmp(a,1) ' / ' SlaveName{2} tmp(a,2) ': '...
                    num2str(meanPct(a)*100, '%.1f') '±' num2str(stdPct(a)*100, '%.1f') '%'];
            end
    end
end
if ip.Results.ShowLegend
    fpos = [2 2 17 5.5];
    aposy = 1.5;
else
    aw = ceil(na/ah);
    fpos = [2 2 8+7*(aw-1) 6+4.5*(ah-1)];
    aposy = 2;
end

A = cell(nCh,nc);

ha = zeros(na,1);
figure(fset.fOpts{:}, 'Position', fpos, 'Name', 'Intensity cohorts');
for a = 1:na
    y0 = ceil(a/2);
    x0 = 1-mod(a,2);
    ha(a) = axes(fset.axOpts{:}, 'Position', [1.5+x0*6.75 aposy+(ah-y0)*4.25 6 3.5]); %#ok<LAXES>
    hold on;
    
    % combination for these axes: sigCombIdx(a)
    for i = 1:nd
        for c = 1:nc
            switch nCh
                case 1
                    res(i).sigComb{c} = res(i).sigIdx{c}==1;
                case 2
                    res(i).sigComb{c} = res(i).sigIdx{c}(:,2)==sigCombIdx(a,1);
                case 3
                    res(i).sigComb{c} = res(i).sigIdx{c}(:,2)==sigCombIdx(a,1) &...
                        res(i).sigIdx{c}(:,3)==sigCombIdx(a,2);
            end
        end
    end
    
    for c = nc:-1:1
        for ch = chVec; % plot master channel last
            if nd > 1
                % means for each data set
                AMat = arrayfun(@(x) ip.Results.AvgFun(x.interpTracks{ch,c}(x.sigComb{c},:),1), res, 'UniformOutput', false);
                AMat = vertcat(AMat{:});
                A{ch,c} = nanmean(AMat,1);
                SEM = nanstd(AMat,[],1)/sqrt(nd);
                Amin = A{ch,c} - SEM;
                Aplus = A{ch,c} + SEM;
            else
                % if input is a single data set, show median + percentiles
                M = prctile(res(1).interpTracks{ch,c}(res(1).sigComb{c},:), [25 50 75], 1);
                A{ch,c} = M(2,:);
                Amin = M(1,:);
                Aplus = M(3,:);
            end
            if ip.Results.ShowVariation
                fill([cT{c} cT{c}(end:-1:1)], sf(ch)*[Amin Aplus(end:-1:1)], cv{ch}(c,:), 'EdgeColor', cmap{ch}(c,:));
            end
            if ~ip.Results.FrontLayer
                plot(cT{c}, sf(ch)*A{ch,c}, ip.Results.LineStyle, 'Color', cmap{ch}(c,:), 'LineWidth', 1);
            end
            cohorts.t{c} = cT{c};
            cohorts.Amin{ch,c} = Amin;
            cohorts.Aplus{ch,c} = Aplus;
            cohorts.A{ch,c} = A{ch,c};
        end
    end
    
    for ch = chVec
        % Plot mean/median in front
        if ip.Results.FrontLayer
            for c = nc:-1:1
                plot(cT{c}, sf(ch)*A{ch,c}, ip.Results.LineStyle, 'Color', cmap{ch}(c,:), 'LineWidth', 1);
            end
        end
        
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
    
end

if isempty(YLim)
    YLim = get(ha, 'YLim');
    if na>1
        YLim = vertcat(YLim{:});
        YLim = [min(YLim(:,1)) max(YLim(:,2))];
    end
end

if ~isempty(sigCombIdx)
    for a = 1:na
        text(XLim(2), YLim(2), atext{a}, fset.sfont{:},...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'Parent', ha(a));
    end
end

if ip.Results.ShowPct && nCh>2
    pos = get(gcf, 'Position');
    pos(3) = pos(3)+3;
    set(gcf, 'Position', pos);
    dy = 0.75;
    hav = axes(fset.axOpts{:}, 'Position', [1.5+ceil(na/ah)*6.75 aposy+ah*3.5+(ah-1)*dy-1.6 2.4 1.6]);
    vennplot(meanPct(2), meanPct(3), meanPct(1), ip.Results.SlaveName,...
        'Handle', hav, 'Font', fset.sfont);
    axis off;
end
    
   
%     if ip.Results.ShowLegend
%         cohortLabels = arrayfun(@(i) [' ' num2str(cohortBounds(i)) '-' num2str(cohortBounds(i+1)-framerate) ' s'], 1:nc, 'UniformOutput', false);
%         hl = legend(hp, [cohortLabels cohortLabels], 'Location', 'SouthEast');
%         set(hl, 'Box', 'off', fset.tfont{:}, 'Position', [6.75+7.65 1.5 1.25 3.5]);
%     end
    
%     if ip.Results.ShowPct
%         axes(fset.axOpts{:}, 'Position', [15.5 2 3 2.5], 'TickLength', fset.TickLength*6/3);
%         barplot2(mean(M,1)', std(M,[],1)', 'Angle', 0, 'BarWidth', 1, 'GroupDistance', 1,...
%             'FaceColor', 0.8*[1 1 1], 'EdgeColor', 0.4*[1 1 1],...
%             'YLim', [0 100], 'LineWidth', 1);
%         set(gca, 'FontSize', 8);
%         
%         h = title(['% ' ip.Results.SlaveName ' pos. CCPs'], fset.sfont{:});
%         %h = ylabel('% CCPs/cohort', fset.lfont{:});
%         pos = get(h, 'Position');
%         %pos(1) = 0.8*pos(1);
%         pos(2) = 1.1*pos(2);
%         set(h, 'Position', pos);
%         set(gca, 'YTick', 0:20:100, 'XTickLabel', cohortLabels);
%         rotateXTickLabels(gca, 'AdjustFigure', false);
%         xlabel('Lifetime cohort', fset.lfont{:});
%     end
    idx = 1-mod(1:na,2)~=0;
    set(ha(idx), 'YTickLabel', []);
    arrayfun(@(x) ylabel(x, 'Fluo. intensity (A.U.)', fset.lfont{:}), ha(~idx));
% end

set(ha, 'XLim', XLim, 'XTick', XTick, 'YLim', YLim);
if ~isempty(ip.Results.YTick)
    set(ha, 'YTick', ip.Results.YTick);
end

if ip.Results.ShowLegend
    arrayfun(@(x) xlabel(x, 'Time (s)', fset.lfont{:}), ha);
else
    idx = ah-ceil((1:na)/2)>0;
    set(ha, 'XTick', XTick, 'XTickLabel', cohortLabels);
    set(ha(idx), 'XTickLabel', []);
    for i = find(~idx)
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


