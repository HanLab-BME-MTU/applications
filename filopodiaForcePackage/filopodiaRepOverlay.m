function repInfo = filopodiaRepOverlay(MLAll, nameList, figPath, varargin)
%FILOPODIAREPOVERLAY  Pick a representative movie per condition and draw the
%filopodia segmentation (body outline + tip points + tip->base shaft lines)
%on top of the talin channel image, for visual QC of detection quality.
%
%   repInfo = filopodiaRepOverlay(MLAll, nameList, figPath)
%   repInfo = filopodiaRepOverlay(MLAll, nameList, figPath, 'ChannelIndex',1, ...
%                'Frame','mid', 'Select','medianCount')
%
% INPUTS
%   MLAll    : array of MovieList (one per condition), as loaded in filopodiaForceBatch
%   nameList : cellstr of condition names (same order as MLAll)
%   figPath  : folder to save the overlay figure
%
% PARAMETERS
%   'ChannelIndex' (1)        talin channel
%   'Frame'        ('mid')    'mid' | 'first' | integer frame index
%   'Select'       ('medianCount')  'medianCount' = movie whose filopodia
%                              count is closest to the condition median
%                              (most representative); 'maxCount' = busiest movie
%
% OUTPUT
%   repInfo(ii).condition / .moviePath / .movieIndex / .frame / .nFilo
%              prints the MD path so you can open it in movieViewer if wanted.
% Sangyoon J. Han / 2026

ip = inputParser;
ip.addParameter('ChannelIndex',1,@isnumeric);
ip.addParameter('Frame','mid');
ip.addParameter('Select','medianCount',@ischar);
ip.parse(varargin{:});
iChan   = ip.Results.ChannelIndex;
frameSel= ip.Results.Frame;
selMode = ip.Results.Select;

nc = numel(MLAll);
repInfo = struct('condition',{},'moviePath',{},'movieIndex',{},'frame',{},'nFilo',{});

hFig = figure('Color','w','Position',[100 100 460*nc 460]);

for ii = 1:nc
    movies = MLAll(ii).movies_;
    nMov = numel(movies);

    % --- count filopodia per movie (from P4) to pick representative ---
    cnt = nan(nMov,1);
    filoCell = cell(nMov,1);
    for k = 1:nMov
        [filo,~] = loadP4(movies{k}, iChan);
        if isempty(filo), continue; end
        filoCell{k} = filo;
        cnt(k) = numel(filo);
    end
    valid = find(~isnan(cnt));
    if isempty(valid)
        warning('Condition %s: no P4 output in any movie.', nameList{ii});
        subplot(1,nc,ii); axis off; title(sprintf('%s: no data',nameList{ii}));
        continue;
    end
    switch selMode
        case 'maxCount'
            [~,rel] = max(cnt(valid));
        otherwise   % medianCount
            med = median(cnt(valid));
            [~,rel] = min(abs(cnt(valid)-med));
    end
    kSel = valid(rel);
    MD   = movies{kSel};
    filo = filoCell{kSel};

    % --- choose frame ---
    if ischar(frameSel)
        if strcmp(frameSel,'first'), fr = 1; else, fr = max(1,round(MD.nFrames_/2)); end
    else
        fr = min(max(1,round(frameSel)),MD.nFrames_);
    end

    % --- load talin image + body mask for this frame ---
    img = double(MD.channels_(iChan).loadImage(fr));
    bodyMask = loadBodyMask(MD, fr);

    % --- draw ---
    ax = subplot(1,nc,ii);
    imshow(img,[],'Parent',ax); hold(ax,'on');
    colormap(ax,gray);
    % body outline
    if ~isempty(bodyMask)
        B = bwboundaries(bodyMask);
        for b = 1:numel(B)
            plot(ax, B{b}(:,2), B{b}(:,1), 'c-', 'LineWidth',1);
        end
    end
    % filopodia: tip point + tip->base line, for filopodia present in this frame
    nShown = 0;
    for f = 1:numel(filo)
        j = find(filo(f).frames==fr,1);
        if isempty(j), continue; end
        tp = filo(f).tipPos(j,:); bp = filo(f).basePos(j,:);
        plot(ax,[bp(1) tp(1)],[bp(2) tp(2)],'-','Color',[1 0.55 0],'LineWidth',1.2);
        plot(ax,tp(1),tp(2),'.','Color',[1 1 0],'MarkerSize',9);
        plot(ax,bp(1),bp(2),'.','Color',[0 1 0],'MarkerSize',7);
        nShown = nShown + 1;
    end
    hold(ax,'off');
    title(ax,sprintf('%s  (movie %d, frame %d)\n%d filo total, %d in frame', ...
        nameList{ii}, kSel, fr, numel(filo), nShown), 'Interpreter','none','FontSize',9);

    % --- record + print path for movieViewer ---
    mp = '';
    try, mp = MD.getFullPath; catch, try, mp = fullfile(MD.movieDataPath_,MD.movieDataFileName_); catch, end, end
    repInfo(ii).condition  = nameList{ii};
    repInfo(ii).moviePath  = mp;
    repInfo(ii).movieIndex = kSel;
    repInfo(ii).frame      = fr;
    repInfo(ii).nFilo      = numel(filo);
    fprintf('[%s] representative = movie %d (count=%d, median basis). Frame %d.\n', ...
        nameList{ii}, kSel, numel(filo), fr);
    fprintf('     movieViewer path: %s\n', mp);
end

% legend-ish annotation
annotation(hFig,'textbox',[0.0 0.0 1 0.05],'String', ...
    'cyan = cell body outline   |   orange = tip\rightarrowbase shaft   |   yellow dot = tip   |   green dot = base', ...
    'EdgeColor','none','HorizontalAlignment','center','FontSize',9);

% --- save ---
if exist(figPath,'dir')~=7, mkdir(figPath); end
fbase = fullfile(figPath,'representativeSegmentation');
try
    hgexport(hFig, fbase, hgexport('factorystyle'), 'Format','eps');
    hgsave(hFig, fbase, '-v7.3');
catch
end
print(hFig, [fbase '.tif'], '-dtiff', '-r150');
fprintf('Saved representative overlay -> %s.tif\n', fbase);
end

% =====================================================================
function [filo, MD] = loadP4(MD, iChan)
filo = [];
try
    ip = MD.getProcessIndex('FilopodiaClassificationProcess',1,0);
    if isempty(ip), return; end
    pr = MD.processes_{ip};
    f = pr.outFilePaths_{1,iChan};
    if isempty(f)||exist(f,'file')~=2
        f = fullfile(pr.funParams_.OutputDirectory,'filoClassification.mat');
    end
    if exist(f,'file')~=2, return; end
    S = load(f,'filopodia');
    filo = S.filopodia;
catch
end
end

% =====================================================================
function bodyMask = loadBodyMask(MD, fr)
bodyMask = [];
try
    ip = MD.getProcessIndex('FilopodiaSegmentationProcess',1,0);
    if isempty(ip), return; end
    pr = MD.processes_{ip};
    od = pr.outFilePaths_{1, find(~cellfun(@isempty,pr.outFilePaths_(1,:)),1)};
    if isempty(od) || exist(od,'dir')~=7
        od = pr.funParams_.OutputDirectory;
    end
    fn = fullfile(od, sprintf('filoSeg_frame_%04d.mat', fr));
    if exist(fn,'file')~=2, return; end
    S = load(fn,'bodyMask');
    bodyMask = S.bodyMask;
catch
end
end