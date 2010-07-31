function makeTropoFigure4(paths, outputDirectory)

colors = [
   0.983333333333333   1.000000000000000   0.800000000000000;
   0.360000000000000   0.630000000000000   0.900000000000000;
   0.060000000000000   0.330000000000000   0.600000000000000;
   0.700000000000000   0.245000000000000   0.245000000000000;
   0.550000000000000                   0                   0;
   0.250000000000000                   0                   0; ];

nMovies = numel(paths);

% Panel A
dataA1 = cell(nMovies,1);
dataA2 = cell(nMovies,1);

% Panel B
dataB = cell(nMovies,2);

for iMovie = 1:nMovies
    % Load Movie data
    load(fullfile(paths{iMovie}, 'movieData.mat'));
    
    pixelSize = movieData.pixelSize_nm;
    %timeInterval = movieData.timeInterval_s;

    % Load density scores
    load(fullfile(movieData.density.directory, movieData.density.channelDirectory{1}, 'densityScores.mat'));
    
    distFromEdge = vertcat(densityScores(:).distFromEdge) * pixelSize; %#ok<NODEF>
    averageDensity = vertcat(densityScores(:).averageDensity) * pixelSize;
    protrusionState = vertcat(densityScores(:).protrusionState);
    protrusionPersistence = vertcat(densityScores(:).protrusionPersistence);
    %protrusionSpeed = abs(vertcat(densityScores(:).protrusionSpeed) * pixelSize / timeInterval);
    
    %
    % Data for Panel A1
    %
    
    maxDistFromEdge = min(15000,max(distFromEdge));    
    dist = 0:500:maxDistFromEdge;
    
    dataA1{iMovie} = zeros(numel(dist)-1,1);

    for i = 1:numel(dist)-1
        dataA1{iMovie}(i) = ...
            mean(averageDensity(distFromEdge > dist(i) & distFromEdge <= dist(i+1) & protrusionState == 2));
    end
    
    %
    % Data for Panel A2
    %

    dataA2{iMovie} = zeros(numel(dist)-1,1);
    
    for i = 1:numel(dist)-1
        dataA2{iMovie}(i) = ...
            mean(averageDensity(distFromEdge > dist(i) & distFromEdge <= dist(i+1) & protrusionState == 3));
    end
    
    %
    % Data for Panel B
    %
    
%     maxProtPersistence = max(protrusionPersistence(protrusionState == 2));
%     maxRetPersistence = max(protrusionPersistence(protrusionState == 3));
%     
%     protPersistence = 1:maxProtPersistence;
%     retPersistence = 1:maxRetPersistence;
%     
%     dataA1Lifetime{iMovie} = zeros(numel(protPersistence), 1);
%     
%     for i = 1:numel(protPersistence)
%         dataA1Lifetime{iMovie}(i) = ...
%             mean(averageDensity(distFromEdge < 2000 & protrusionPersistence == protPersistence(i)));
%     end

    dataB{iMovie,1} = ...
        averageDensity(distFromEdge < 2000 & protrusionState == 2 & protrusionPersistence == 1);

    dataB{iMovie,2} = ...
        averageDensity(distFromEdge < 2000 & protrusionState == 3 & protrusionPersistence == 1);
    
    %
    % Data for Panel B2
    %
    
%     dataA2Lifetime{iMovie} = zeros(numel(retPersistence), 1);
%     
%     for i = 1:numel(retPersistence)
%         dataA2Lifetime{iMovie}(i) = ...
%             mean(averageDensity(distFromEdge < 2000 & protrusionPersistence == retPersistence(i)));
%     end
    
end

%
% Panel A1
%

% convert dist in microns:
dist = dist / 1000;

hFig = figure('Visible', 'off');
set(gca,'XTick',dist(1:4:end-1));
set(gca, 'FontName', 'Helvetica', 'FontSize', 18);
set(gcf, 'Position', [680 678 650 450], 'PaperPositionMode', 'auto');

averageDensityTotal = horzcat(dataA1{:});

% plot x axis in um
h = line(dist(1:end-1), averageDensityTotal,'LineWidth',2);

for i=1:numel(h)
    set(h(i),'Color', colors(i,:));
end

h = legend({'TM2', 'TM4', 'TM5'}); legend('boxoff');
hC = get(h,'Children');
hC = hC(1:2:end);

for i=1:numel(hC)
    hCC = get(hC(i),'Children');
    set(hCC,'FaceColor', colors(numel(hC) - i + 1,:));
end

xlabel(['Distance away from cell edge during protrusion (' char(181) 'm)']);
ylabel('Speckle Density (nm^{-2})');

fileName = [outputDirectory filesep 'Fig4_A1.eps'];
print(hFig, '-depsc', fileName);
fixEpsFile(fileName);
close(hFig);

%
% Panel A2
%

hFig = figure('Visible', 'off');
set(gca,'XTick',dist(1:4:end-1));
set(gca, 'FontName', 'Helvetica', 'FontSize', 18);
set(gcf, 'Position', [680 678 650 450], 'PaperPositionMode', 'auto');

averageDensityTotal = horzcat(dataA2{:});

% plot x axis in um
h = line(dist(1:end-1), averageDensityTotal,'LineWidth',2);

for i=1:numel(h)
    set(h(i),'Color', colors(i,:));
end

h = legend({'TM2', 'TM4', 'TM5'}); legend('boxoff');
hC = get(h,'Children');
hC = hC(1:2:end);

for i=1:numel(hC)
    hCC = get(hC(i),'Children');
    set(hCC,'FaceColor', colors(numel(hC) - i + 1,:));
end

xlabel(['Distance away from cell edge during retraction (' char(181) 'm)']);
ylabel('Speckle Density (nm^{-2})');

fileName = [outputDirectory filesep 'Fig4_A2.eps'];
print(hFig, '-depsc', fileName);
fixEpsFile(fileName);
close(hFig);

%
% Panel B
%

hFig = figure('Visible', 'off');
set(gca, 'FontName', 'Helvetica', 'FontSize', 18);
set(gcf, 'Position', [680 678 650 450], 'PaperPositionMode', 'auto');

mu = cellfun(@mean, dataB);
sigma = cellfun(@std, dataB);
h = bar(gca, mu', 'group'); hold on;

XTicks = zeros(size(mu));

for i=1:numel(h)
    hC = get(h(i), 'Children');
    set(hC,'FaceColor', colors(i,:));
    XData = get(hC, 'XData');
    XTicks(i,1:2) = .5 * (XData(1,1:2) + XData(3,1:2));
end

errorbar(XTicks(:),mu(:),sigma(:),'xk'); hold off;

set(gca,'XTickLabel',{'Protrusion', 'Retraction'})
ylabel('D_0 (nm^{-2})');

h = legend({'TM2', 'TM4', 'TM5NM1'}); legend('boxoff');
hC = get(h,'Children');
hC = hC(1:2:end);
for i=1:numel(hC)
    hCC = get(hC(i),'Children');
    set(hCC,'FaceColor', colors(numel(hC) - i + 1,:));
end

fileName = [outputDirectory filesep 'Fig4_B.eps'];
print(hFig, '-depsc', fileName);
fixEpsFile(fileName);
close(hFig);

% %
% % Panel B2
% %
% 
% hFig = figure('Visible', 'off');
% set(gca, 'FontName', 'Helvetica', 'FontSize', 18);
% set(gcf, 'Position', [680 678 650 450], 'PaperPositionMode', 'auto');
% 
% % append with zero
% maxLength = max(cellfun(@numel, dataA2Lifetime));
% dataA2Lifetime = cellfun(@(x)...
%     padarray(x, maxLength - numel(x), 0, 'post'), ...
%     dataA2Lifetime, 'UniformOutput', false);
% 
% averageDensityTotal = horzcat(dataA2Lifetime{:});
% 
% % plot x axis in um
% h = bar(1:size(averageDensityTotal,1), averageDensityTotal, 'group');
% 
% for i=1:numel(h)
%     hC = get(h(i), 'Children');
%     set(hC,'FaceColor', colors(i,:));
% end
% 
% h = legend({'TM2', 'TM4', 'TM5'}); legend('boxoff');
% hC = get(h,'Children');
% hC = hC(1:2:end);
% 
% for i=1:numel(hC)
%     hCC = get(hC(i),'Children');
%     set(hCC,'FaceColor', colors(numel(hC) - i + 1,:));
% end
% 
% xlabel('Retraction Lifetime (#frame)');
% ylabel('Speckle Density (nm^{-2})');
% 
% fileName = [outputDirectory filesep 'Fig4_B2.eps'];
% print(hFig, '-depsc', fileName);
% fixEpsFile(fileName);
% close(hFig);

