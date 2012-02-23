function plotInitIntensityVsLifetime(data, T)


nd = numel(data);

for k = 1:nd
    
    load([data(k).source 'Tracking' filesep 'trackAnalysis.mat']);
    [tracksAcc tracksRej] = runTrackPostProcessing(data(k), tracks);
    
    % tracks used for threshold
    tTracks = tracksAcc([tracksAcc.lifetime_s]==T);
    
    maxA = arrayfun(@(i) max(i.A), tTracks);
    meanMaxA = mean(maxA);
    
    selTracks = tracksAcc([tracksAcc.lifetime_s]>T+1);
    
    % lifetime (in frames) until A > T
    idelta = arrayfun(@(i) find(i.A>meanMaxA, 1, 'first'), selTracks, 'UniformOutput', false);
    maxA = arrayfun(@(i) max(i.A), selTracks);
    
    valid = cellfun(@(i) ~isempty(i), idelta);
    idelta = idelta(valid);
    imax = maxA(valid);
    ilft = [selTracks(valid).lifetime_s];
    
    %idelta(didx==0) = {NaN};
    %idelta = [idelta{:}];
    delta{k} = [idelta{:}];
    maxInt{k} = imax;
    lft{k} = ilft;
    
    %---------------------------------
    lv = 4:40;
    for i = 1:numel(lv)
        itracks = tracksAcc([tracksAcc.lifetime_s]==lv(i));
        maxA = arrayfun(@(i) max(i.A), itracks);
        meanMaxA = mean(maxA);
        idelta = arrayfun(@(i) find(i.A>meanMaxA, 1, 'first'), tracksAcc([tracksAcc.lifetime_s]>lv(k)+1), 'UniformOutput', false);
        valid = cellfun(@(i) ~isempty(i), idelta);
        idelta = idelta(valid);
        idelta = [idelta{:}];
        
        ni = hist(idelta,1:100);
        dfVect{k}(i) = find(ni==max(ni), 1, 'first');
        %dfVect(k) = find(ni==round(nanmean(delta)), 1, 'first');
        
    end
    
end
% pool all data sets
delta = [delta{:}];
maxInt = [maxInt{:}];
lft = [lft{:}];

cutoff_f = 4;
N = 103;
w = N./(N-cutoff_f+1:-1:1);
ni = hist(delta, 1:100);
ni = ni.*w;
ni = ni/sum(ni);
%%

dfMean = mean(vertcat(dfVect{:}),1);
dfStd = std(vertcat(dfVect{:}),[],1);


figure; 
hold on;
% plot(0:40,0:40,'r');
plot(lv, dfMean, 'k.-');
plot(lv, dfMean+dfStd, 'r-');
plot(lv, dfMean-dfStd, 'r-');

axis equal tight
axis([1 lv(end) 1 max(dfMean)+1]);

% x frames needed to surpass max intensity of length==x tracks


%%
pos = get(0, 'DefaultFigurePosition');
pos(4) = 600;
figure('Position', pos, 'PaperPositionMode', 'auto', 'Units', 'pixels', 'Color', [1 1 1]); 
xi = 1:100;
%0.1300    0.1100    0.7750    0.8150

axes('Position', [0.13 0.7+0.06 0.7750 0.2]);
hold on;
set(gca, 'XTick', []);
intDist = arrayfun(@(i) maxInt(delta==i), 1:100, 'UniformOutput', false);
meanInt = cellfun(@(i) mean(i), intDist);
stdInt = cellfun(@(i) std(i), intDist);

plot(xi, meanInt);
plot(xi, meanInt-stdInt);
plot(xi, meanInt+stdInt);
set(gca, 'XLim', [0 40]);
ylabel('Max. intensity');

axes('Position', [0.13 0.5+0.03 0.7750 0.2])
% set(gca, 'XTick', []);
% set(gca, 'TickLength', [0 0], 'XTickLabel', []);
hold on;
% grid on;
lftDist = arrayfun(@(i) lft(delta==i), 1:100, 'UniformOutput', false);
medLft = cellfun(@(i) median(i), lftDist);
% stdLft = cellfun(@(i) std(i), lftDist);

p10Lft = zeros(1,100);
p90Lft = zeros(1,100);
for i = 1:100
    if numel(lftDist{i})>5
        [f_ecdf, i_ecdf] = ecdf(lftDist{i});
        f_ecdf = f_ecdf(2:end);
        i_ecdf = i_ecdf(2:end);
        ival = interp1(f_ecdf, i_ecdf, [0.1 0.9]);
    else
        ival = NaN(1,2);
    end
    p10Lft(i) = ival(1);
    p90Lft(i) = ival(2);
end

plot(xi, medLft);
plot(xi, p10Lft);
plot(xi, p90Lft);
set(gca, 'XLim', [0 40]);
ylabel('Lifetime (s)');

axes('Position', [0.13 0.11 0.7750 0.39])

% hist(delta,1:100);
bar(xi, ni, 'BarWidth', 1);
box off;

set(gca, 'XLim', [0 40]);
xlabel(['Frames until ''A'' > mean(max(A_{' num2str(T) '}))'], 'FontSize', 20);
set(gca, 'LineWidth', 2, 'TickDir', 'out', 'Layer', 'top', 'FontSize', 18, 'XTick', 0:5:50);
%%


