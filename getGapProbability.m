% Francois Aguet, 06/13/2012

function gapPDF = getGapProbability(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addOptional('ch', 1);
ip.addParamValue('CohortBounds_s', [10 20 40 60 80 100 120]);
ip.addParamValue('MaxIntensityThreshold', 0);
ip.addParamValue('Display', 'on', @(x) any(strcmpi(x, {'on', 'off'})));
ip.parse(data, varargin{:});
cohortBounds = ip.Results.CohortBounds_s;

nc = numel(cohortBounds)-1;

lftData = getLifetimeData(data);
nd = numel(data);
framerate = data(1).framerate;
mCh = strcmp(data(1).source, data(1).channels);

% since gaps are SNR-dependent, process data sets separately


% # data points in cohort (including buffer frames)
% iLength = arrayfun(@(c) floor(mean(cohortBounds([c c+1]))/framerate), 1:nc);
iLength = arrayfun(@(c) cohortBounds(c+1)/framerate, 1:nc);
% time vectors for cohorts
cT = arrayfun(@(i) (0:i-1)*framerate, iLength, 'UniformOutput', false);
res = struct([]);
for i = 1:nd
    lifetime_s = lftData(i).lifetime_s([lftData(i).catIdx]==1);
    trackLengths = lftData(i).trackLengths([lftData(i).catIdx]==1);
    maxA = max(lftData(i).intMat_Ia(:,:,mCh), [], 2)';
    
    for c = 1:nc
        cidx = find(cohortBounds(c)<=lifetime_s & lifetime_s<cohortBounds(c+1) & maxA > ip.Results.MaxIntensityThreshold);
        nt = numel(cidx);
        cLengths = trackLengths(cidx);
        interpGaps = zeros(nt,iLength(c));

        % loop through tracks within cohort
        for t = 1:nt
            gaps = double(lftData(i).gapMat_Ia(cidx(t),1:cLengths(t)));
            xi = linspace(1,cLengths(t), iLength(c));
            interpGaps(t,:) = binterp(gaps, xi);
        end
        interpGaps = sum(interpGaps,1);
        %interpGaps = interpGaps/sum(interpGaps);
        interpGaps = interpGaps/nt; % normalize by # tracks
        res(i).gapPDF{c} = interpGaps;
    end
end

cmap = jet(nc);
figure;
hold on;
gapPDF = cell(1,nc);
for c = 1:nc
    % concatenate and average
    G = arrayfun(@(i) i.gapPDF{c}, res, 'UniformOutput', false);
    G = vertcat(G{:});
    mu = mean(G,1);
    mu = mu/sum(mu);
    mu(mu<0) = 0; % correct for numerical errors
    plot(cT{c}, mu, 'Color', cmap(c,:));
    gapPDF{c} = mu;
end
