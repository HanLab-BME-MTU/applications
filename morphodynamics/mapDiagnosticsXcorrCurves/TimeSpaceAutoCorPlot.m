function [acmapOf_, corrMatOf_, meanACOf_] ...
         = TimeSpaceAutoCorPlot(map, mapName, MDtimeInterval_, varargin)
% TimeSpaceAutoCorPlot Compute/draw auto-correlation map, correlation
% matrix map, and their means over space and time separately.
% Jungsik Noh, 2016/10/05

ip = inputParser;
ip.addParameter('figFlag', 'off');
parse(ip, varargin{:})
p = ip.Results;
figFlag = p.figFlag;


%%
 
xx1 = 1:ceil(size(map, 2)/2);
xx2 = 1:ceil(size(map, 1)/2);


%%%% over time
acmapThresh = 2/sqrt(size(map, 2));
acmap = autoCorrMap(map, 'maxLag', xx1(end));  
 
%
corAvg_autocor = mean(acmap(:, xx1), 1, 'omitnan');

%%%%
%%%% over space
spatialAcmapThresh = 2/sqrt(size(map, 1));
spatialAcmap = autoCorrMap(map', 'maxLag', xx2(end));
%
corAvg_spatialAutocor = mean(spatialAcmap(:, xx2), 1, 'omitnan');

%
acmap1 = acmap(:, 2:size(acmap, 2));
spatialAcmap1 = spatialAcmap(:, 2:size(spatialAcmap, 2));

%%%% subplots

acmapOf_ = figure('Visible', figFlag);
subplot(2, 1, 1);
figtmp = imagesc(acmap1, [-1 1]);
colorbar;colormap(jet)
figtmp.AlphaData = 1-isnan(acmap1);

axis xy;xlabel('Time lag (s)');ylabel('Window')
ax = gca;
curTick = ax.XTick;
ax.XTickMode = 'manual';
ax.XTick = curTick;
ax.XTickLabel = curTick*MDtimeInterval_;

set(gca, 'XGrid', 'on')
title(['AutoCorr of ', mapName, ' (Thresh: ', num2str(acmapThresh, 4), ')'])

%
subplot(2, 1, 2);
figtmp = imagesc(spatialAcmap1, [-1 1]);
colorbar;colormap(jet)
figtmp.AlphaData = 1-isnan(spatialAcmap1);

axis xy;xlabel('Window lag');ylabel('Time (s)')
ax = gca;
curTick = ax.YTick;
ax.YTickMode = 'manual';
ax.YTick = curTick;
ax.YTickLabel = (curTick-1)*MDtimeInterval_;

set(gca, 'XGrid', 'on')
title(['spatial AutoCorr of ', mapName, ' (Thresh: ', num2str(spatialAcmapThresh, 4), ')'])

%



%  corr(x_t1, x_t2) along windows
%   corr(x^w1, x^w2) along time frames

% corr(x_t1, x_t2) along windows
% Deal with nan's
cormatTime = corr(map, 'rows', 'pairwise');                     

xx3 = 0:(length(corAvg_autocor)-1);

corAvg3 = [];
for k = xx3
    tmp = mean(diag(cormatTime, k), 'omitnan');
    corAvg3 = [corAvg3 tmp];
end
%corAvg3;

% corr(x^w1, x^w2) along time frames
% Deal with nan's
cormatWin = corr(map', 'rows', 'pairwise');

xx4 = 0:(length(corAvg_spatialAutocor)-1);

corAvg4 = [];
for k = xx4
    tmp = mean(diag(cormatWin, k), 'omitnan');
    corAvg4 = [corAvg4 tmp];
end
%corAvg4;


corrMatOf_ = figure('Visible', figFlag);                             % Figure name
subplot(1,2,1);
figtmp = imagesc(cormatTime, [-1 1]);
colorbar;colormap(jet)
figtmp.AlphaData = 1-isnan(cormatTime);
title({mapName, 'corr(x_{t1}, x_{t2})'})
xlabel('Time (s)');ylabel('Time (s)')

ax = gca;
curTick = ax.XTick;
ax.XTickMode = 'manual';
ax.XTick = curTick;
ax.XTickLabel = (curTick-1)*MDtimeInterval_;

curTick = ax.YTick;
ax.YTickMode = 'manual';
ax.YTick = curTick;
ax.YTickLabel = (curTick-1)*MDtimeInterval_;



subplot(1,2,2)
figtmp = imagesc(cormatWin, [-1 1]);
colorbar;colormap(jet)
figtmp.AlphaData = 1-isnan(cormatWin);

title({mapName, 'corr(x^{w1}, x^{w2})'})
xlabel('Window');ylabel('Window')


%%%%    the cases of 'along with windows' tend to be larger 
%%%%    => (maybe due to) window wise heterogeniety is larger

% mean of auto correlations from 2 approaches
tlag = xx3*MDtimeInterval_;

meanACOf_ = figure('Visible', figFlag);                                  % Figure name
subplot(2,1,1)
plot(tlag, corAvg_autocor, tlag, corAvg3, '--o')
xlabel('Time lag (s)')
ylabel('Avg. auto correlation')
legend('from corrMap', 'from corrMatrix')   %, 'Location','northoutside','Orientation','horizontal')
title(mapName)
h = refline(0, 0);
h.Color = [.5 .5 .5];

subplot(2,1,2)
plot(xx4, corAvg_spatialAutocor, xx4, corAvg4, '--o')
xlabel('Window lag')
ylabel('Avg. auto correlation')
legend('from corrMap', 'from corrMatrix')   %, 'Location','northoutside','Orientation','horizontal')
title(mapName)
h = refline(0, 0);
h.Color = [.5 .5 .5];


end


%%  EOF
%%%%%%%%%
