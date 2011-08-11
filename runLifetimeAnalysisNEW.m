function data = runLifetimeAnalysisNEW(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParamValue('FileName', 'trackAnalysis.mat', @ischar);
ip.addParamValue('Tracks', []);

ip.parse(data, varargin{:});


for k = 1:length(data)
    if isempty(ip.Results.Tracks)
        ta = load([data(k).source 'Tracking' filesep ip.Results.FileName]);
        tracks = ta.tracks;
        data(k).tracks = tracks([tracks.valid]==1);
    else
        data(k).tracks = ip.Results.Tracks;
    end        
    data(k).lifetimes = [data(k).tracks.lifetime_s];
end

% maxLifetime = max([data.lifetimes]);


% Extend all to max. movie length, in case of mismatch
Nmax = max([data.movieLength])-2;

% generate lifetime histograms
for k = 1:length(data)
    dt = data(k).framerate;
    N = data(k).movieLength-2;
    t = (1:N)*dt;
    lftHist = hist(data(k).lifetimes, t);
    
    % apply correction
    % longest observable lifetime: N = movieLength-2
    % relative probabilities:
    % P(obs. lifetime==1) = N
    % P(obs. lifetime==N) = 1
    % => weighting:
    lftHist = lftHist .* N./(N:-1:1);
    
    % Extend
    if N<Nmax
        lftHist = [lftHist Nmax-N];
    end
    
    % Normalization
    lftHist = lftHist / sum(dt*lftHist);
    
    data(k).lftHist = lftHist;
end

%-------------------------
% Mean histogram
%-------------------------
M = vertcat(data.lftHist);
histMean = mean(M,1);
histSEM = std(M,1) / sqrt(length(data));


tfont = {'FontName', 'Helvetica', 'FontSize', 14};
sfont = {'FontName', 'Helvetica', 'FontSize', 18};
lfont = {'FontName', 'Helvetica', 'FontSize', 22};


figure;
hold on;

t = (1:Nmax)*dt;

fill([t t(end:-1:1)], [histMean-histSEM histMean(end:-1:1)+histSEM(end:-1:1)],...
    [1 1 1]*0.7, 'EdgeColor', 'none');
plot(t, histMean, 'k.-');

axis([0 300 0 0.05]);

set(gca, 'LineWidth', 1.5, sfont{:});

xlabel('Lifetime (s)', lfont{:});
ylabel('Frequency', lfont{:});



% TO DO:
% inset with: # used, # persistent, # MS