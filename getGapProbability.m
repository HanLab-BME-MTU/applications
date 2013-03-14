% Francois Aguet, 06/13/2012

function res = getGapProbability(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addOptional('ch', 1);
ip.addParamValue('MaxIntensityThreshold', 0);
ip.addParamValue('Display', 'on', @(x) any(strcmpi(x, {'on', 'off'})));
ip.addParamValue('WindowWidth', 20)
ip.parse(data, varargin{:});

% pool all track gap maps
w = ip.Results.WindowWidth;
lftData = getLifetimeData(data);
mCh = strcmp(data(1).source, data(1).channels);

% align first w frames
gapM = arrayfun(@(i) i.gapMat_Ia(:,1:w,mCh), lftData, 'UniformOutput', false);
gapM = vertcat(gapM{:});

trackLengths = arrayfun(@(i) i.trackLengths(i.catIdx==1), lftData, 'UniformOutput', false);
trackLengths = vertcat(trackLengths{:});

Pstart = sum(gapM(trackLengths>=40,:),1) / size(gapM(trackLengths>=40,:),1);

endM = arrayfun(@(i) i.gapMat_Ia(:,:,mCh), lftData, 'UniformOutput', false);
endM = vertcat(endM{:});
endM = endM(trackLengths>=40,:);
iLengths = trackLengths(trackLengths>=40);
N = size(endM,1);
for k = 1:N
    endM(k,1:w) = endM(k,iLengths(k)-w+1:iLengths(k));
end
endM = endM(:,1:w);
Pend = sum(endM,1) / N;

res.w = w;
res.Pstart = Pstart;
res.Pend = Pend;


% dx = 10;
% xs = 2:w;
% fs = Pstart(2:end);
% xe = w+dx+ (1:w-1);
% fe = Pend(1:end-1);
% 
% % pps = csaps(xs,fs, 0.2);
% % fsi = interpB3SmoothingSpline1D(xs, fs, 0.5);
% 
% figure;
% hold on;
% plot(xs, fs, '.-');
% plot(xe, fe, '.-');
% % plot(xs, fsi, 'r--');
% % plot(xs, pp, 'r--');
% fnplt(pp, 'r--');
% set(gca, 'XTick', [2 10 20 w+dx+[1 9 19]]);


% 
% nd = numel(data);
% framerate = data(1).framerate;
% 
% 
% % since gaps are SNR-dependent, process data sets separately
% 
% 
% % # data points in cohort (including buffer frames)
% % iLength = arrayfun(@(c) floor(mean(cohortBounds([c c+1]))/framerate), 1:nc);
% iLength = arrayfun(@(c) cohortBounds(c+1)/framerate, 1:nc);
% % time vectors for cohorts
% cT = arrayfun(@(i) (0:i-1)*framerate, iLength, 'UniformOutput', false);
% res = struct([]);
% for i = 1:nd
%     lifetime_s = lftData(i).lifetime_s([lftData(i).catIdx]==1);
%     trackLengths = lftData(i).trackLengths([lftData(i).catIdx]==1);
%     maxA = max(lftData(i).A(:,:,mCh), [], 2)';
%     
%     for c = 1:nc
%         cidx = find(cohortBounds(c)<=lifetime_s & lifetime_s<cohortBounds(c+1) & maxA > ip.Results.MaxIntensityThreshold);
%         nt = numel(cidx);
%         cLengths = trackLengths(cidx);
%         interpGaps = zeros(nt,iLength(c));
% 
%         % loop through tracks within cohort
%         for t = 1:nt
%             gaps = double(lftData(i).gapMat_Ia(cidx(t),1:cLengths(t)));
%             xi = linspace(2,cLengths(t)-1, iLength(c));
%             interpGaps(t,:) = binterp(gaps, xi);
%         end
%         interpGaps = sum(interpGaps,1);
%         %interpGaps = interpGaps/sum(interpGaps);
%         interpGaps = interpGaps/nt; % normalize by # tracks
%         res(i).gapPDF{c} = interpGaps;
%     end
% end
% 
% cmap = jet(nc);
% figure;
% hold on;
% gapPDF = cell(1,nc);
% for c = 1:nc
%     % concatenate and average
%     G = arrayfun(@(i) i.gapPDF{c}, res, 'UniformOutput', false);
%     G = vertcat(G{:});
%     mu = mean(G,1);
%     mu = mu/sum(mu);
%     mu(mu<0) = 0; % correct for numerical errors
%     plot(cT{c}, mu, 'Color', cmap(c,:));
%     gapPDF{c} = mu;
% end
