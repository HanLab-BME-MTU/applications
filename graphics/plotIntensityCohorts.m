% Francois Aguet, 05/30/2012

function plotIntensityCohorts(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addOptional('ch', []);
ip.addParamValue('Overwrite', false, @islogical);
% ip.addParamValue('CohortBounds_s', [10 20 40 60 80 100 125 150]);
ip.addParamValue('CohortBounds_s', [10 20 40 60 80 100 120]);
ip.addParamValue('ShowSEM', true, @islogical);
ip.addParamValue('ScaleChannels', 'end', @(x) isempty(x) || any(strcmpi(x, {'start', 'end'})));
ip.parse(data, varargin{:});
cohortBounds = ip.Results.CohortBounds_s;


% if no specific channel is selected, all channels are shown
ch = ip.Results.ch;
mCh = find(strcmp(data(1).source, data(1).channels));


nd = numel(data);
nCh = numel(data(1).channels);
nc = numel(cohortBounds)-1;
sCh = setdiff(1:nCh,mCh);
b = 5;
% loop through data sets, generate cohorts for each
res(1:nd) = struct('cMean', [], 'cSEM', []);
for i = 1:nd
    lftData = getLifetimeData(data(i), 'Overwrite', ip.Results.Overwrite);
    lifetime_s = lftData.lifetime_s([lftData.catIdx]==1);
    trackLengths = lftData.trackLengths([lftData.catIdx]==1);
    
    if ip.Results.ScaleChannels
        % scaling of slave channels relative to master: median of end buffer, median of max. intensity
        medianEndBuffer = squeeze(mean(mean(lftData.([ip.Results.ScaleChannels 'Buffer_Ia']),2),1));
        medianMax = squeeze(mean(max(lftData.intMat_Ia,[],2),1));
        for c = sCh
            lftData.intMat_Ia(:,:,c) = (lftData.intMat_Ia(:,:,c)-medianEndBuffer(c))/medianMax(c)*medianMax(mCh)+medianEndBuffer(mCh);
            lftData.startBuffer_Ia(:,:,c) = (lftData.startBuffer_Ia(:,:,c)-medianEndBuffer(c))/medianMax(c)*medianMax(mCh)+medianEndBuffer(mCh);
            lftData.endBuffer_Ia(:,:,c) = (lftData.endBuffer_Ia(:,:,c)-medianEndBuffer(c))/medianMax(c)*medianMax(mCh)+medianEndBuffer(mCh);
        end
    end
    
    cT = cell(1,nc);
    res(i).cMean = cell(nCh,nc);
    res(i).cSEM = cell(nCh,nc);
    for ch = 1:nCh % channels
        % interpolate tracks to mean cohort length
        for c = 1:nc % cohorts
            % tracks in current cohort
            cidx = find(cohortBounds(c)<=lifetime_s & lifetime_s<cohortBounds(c+1));
            nt = numel(cidx);
            
            % # data points in cohort (with buffer frames)
            iLength = floor(mean(cohortBounds([c c+1]))/data(i).framerate) + 2*b;
            
            interpTracks = zeros(nt,iLength);
            cLengths = trackLengths(cidx);
            % loop through track lengths within cohort
            for t = 1:nt
                A = [lftData.startBuffer_Ia(cidx(t),:,ch) lftData.intMat_Ia(cidx(t),1:cLengths(t),ch) lftData.endBuffer_Ia(cidx(t),:,ch)];
                interpTracks(t,:) = interp1(1:cLengths(t)+2*b, A,...
                    linspace(1,cLengths(t)+2*b, iLength), 'cubic');
                %interpTracks(t,:) = binterp(A, linspace(1,cLengths(t)+2*b, iLength));
            end
            
            cT{c} = (-b:iLength-b-1)*data(i).framerate;
            res(i).cMean{ch,c} = mean(interpTracks,1);
            res(i).cSEM{ch,c} = std(interpTracks,[],1)/sqrt(nt);
        end
    end
end


% cmap = jet(256);
% w = 256/nc/2;
% cmap = interp1(1:256, cmap, linspace(1+w,256-w,nc));

fset = loadFigureSettings();


figure;
hold on;

if nCh==1
    cmap = jet(nc);
    cv = rgb2hsv(cmap);
    cv(:,2) = 0.2;
    cv = hsv2rgb(cv);
    %cmap = ones(nc,3);
    %cmap(:,1) = (nc:-1:1)/nc;  
    %cmap = hsv2rgb(cmap);
    
    for c = nc:-1:1
        if nd>1
            A = arrayfun(@(x) x.cMean{ch,c}, res, 'UniformOutput', false);
            A = vertcat(A{:});
            SEM = std(A,[],1)/sqrt(nd);
            A = mean(A,1);
        else
            A = res(1).cMean{ch,c};
            SEM = res(1).cSEM{ch,c};
        end
        if ip.Results.ShowSEM
            fill([cT{c} cT{c}(end:-1:1)], [A-SEM A(end:-1:1)+SEM(end:-1:1)], cv(c,:), 'EdgeColor', cmap(c,:));
        end
        plot(cT{c}, A, 'Color', cmap(c,:), 'LineWidth', 1.5);
    end
else
    hues = getFluorophoreHues(data(1).markers);
    
    for ch = nCh:-1:1
        trackColor = hsv2rgb([hues(ch) 1 0.8]);
        fillLight = hsv2rgb([hues(ch) 0.4 1]);

        for c = nc:-1:1
            if nd>1
                A = arrayfun(@(x) x.cMean{ch,c}, res, 'UniformOutput', false);
                A = vertcat(A{:});
                SEM = std(A,[],1)/sqrt(nd);
                A = mean(A,1);
            else
                A = res(1).cMean{ch,c};
                SEM = res(1).cSEM{ch,c};
            end
            if ip.Results.ShowSEM
                fill([cT{c} cT{c}(end:-1:1)], [A-SEM A(end:-1:1)+SEM(end:-1:1)], fillLight, 'EdgeColor', trackColor);
            end
            plot(cT{c}, A, 'Color', trackColor, 'LineWidth', 1.5);
        end
    end
end
set(gca, fset.axOpts{:}, 'XLim', [-b*data(1).framerate cohortBounds(end)]);
xlabel('Time (s)', fset.lfont{:});
ylabel('Intensity (A.U.)', fset.lfont{:});
