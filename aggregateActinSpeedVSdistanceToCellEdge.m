function aggregateActinSpeedVSdistanceToCellEdge(conDirectory, creDirectory, outputDirectory, alpha, n)

cmd1 = ['find ' conDirectory ' -name ''*_ActinSpeedOutsideAdhesionVSdistanceToEdge_DATA.mat'''];
cmd2 = ['find ' creDirectory ' -name ''*_ActinSpeedOutsideAdhesionVSdistanceToEdge_DATA.mat'''];
[err1, result1] = system(cmd1);
[err2, result2] = system(cmd2);

assert(~(err1 && err2));

% Read the paths
C = textscan(result1, '%s\n');
conPaths = C{1};
C = textscan(result2, '%s\n');
crePaths = C{1};

% Number of movies
nConMovies = numel(conPaths);
nCreMovies = numel(crePaths);

% Number of pixels in bin 1 and 2 for control and null movies
nConLP = zeros(nConMovies,1);
nConLM = zeros(nConMovies,1);
nCreLP = zeros(nCreMovies,1);
nCreLM = zeros(nCreMovies,1);

% Average Actin speed in bin 1 and 2 for control and null movies
muConLP = zeros(nConMovies,1);
muConLM = zeros(nConMovies,1);
muCreLP = zeros(nCreMovies,1);
muCreLM = zeros(nCreMovies,1);

% Standard deviation of Actin speed in bin 1 and 2 for control and null
% movies
stdConLP = zeros(nConMovies,1);
stdConLM = zeros(nConMovies,1);
stdCreLP = zeros(nCreMovies,1);
stdCreLM = zeros(nCreMovies,1);

totalConLP = cell(nConLP,1);
totalConLM = cell(nConLM,1);
totalCreLP = cell(nCreLP,1);
totalCreLM = cell(nCreLM,1);

for iMovie = 1:nConMovies    
    % load data
    load(conPaths{iMovie});
    
    totalConLP{iMovie} = data{1};
    totalConLM{iMovie} = data{2};
    
    nConLP(iMovie) = numel(data{1}); %#ok<*USENS>
    muConLP(iMovie)= mean(data{1});
    stdConLP(iMovie) = std(data{1});
    
    nConLM(iMovie) = numel(data{2});
    muConLM(iMovie) = mean(data{2});
    stdConLM(iMovie) = std(data{2});
    
    fprintf(1, '%s: (BIN 1) number of pixels = %d, mean = %f, std = %f\n', getDirFromPath(conPaths{iMovie}), nConLP(iMovie), muConLP(iMovie), stdConLP(iMovie));
    fprintf(1, '%s: (BIN 2) number of pixels = %d, mean = %f, std = %f\n', getDirFromPath(conPaths{iMovie}), nConLM(iMovie), muConLM(iMovie), stdConLM(iMovie));    
end

for iMovie = 1:nCreMovies    
    % load data
    load(crePaths{iMovie});
    
    totalCreLP{iMovie} = data{1};
    totalCreLM{iMovie} = data{2};
    
    nCreLP(iMovie) = numel(data{1}); %#ok<*USENS>
    muCreLP(iMovie)= mean(data{1});
    stdCreLP(iMovie) = std(data{1});
    
    nCreLM(iMovie) = numel(data{2});
    muCreLM(iMovie) = mean(data{2});
    stdCreLM(iMovie) = std(data{2});
    
    fprintf(1, '%s: (BIN 1) number of pixels = %d, mean = %f, std = %f\n', getDirFromPath(crePaths{iMovie}), nCreLP(iMovie), muCreLP(iMovie), stdCreLP(iMovie));
    fprintf(1, '%s: (BIN 2) number of pixels = %d, mean = %f, std = %f\n', getDirFromPath(crePaths{iMovie}), nCreLM(iMovie), muCreLM(iMovie), stdCreLM(iMovie));    
end

nTotalConLP = sum(nConLP);
nTotalConLM = sum(nConLM);
nTotalCreLP = sum(nCreLP);
nTotalCreLM = sum(nCreLM);

muTotalConLP = mean(muConLP);
muTotalConLM = mean(muConLM);
stdTotalConLP = sqrt((1/nConMovies) * sum(stdConLP.^2));
stdTotalConLM = sqrt((1/nConMovies) * sum(stdConLM.^2));

muTotalCreLP = mean(muCreLP);
muTotalCreLM = mean(muCreLM);
stdTotalCreLP = sqrt((1/nCreMovies) * sum(stdCreLP.^2));
stdTotalCreLM = sqrt((1/nCreMovies) * sum(stdCreLM.^2));

totalConLP = vertcat(totalConLP{:});
totalConLM = vertcat(totalConLM{:});
totalCreLP = vertcat(totalCreLP{:});
totalCreLM = vertcat(totalCreLM{:});

% sub-sample data:
perms = randperm(numel(totalConLP));
totalConLP = totalConLP(perms(1:min(n,numel(perms))));
perms = randperm(numel(totalConLM));
totalConLM = totalConLM(perms(1:min(n,numel(perms))));
perms = randperm(numel(totalCreLP));
totalCreLP = totalCreLP(perms(1:min(n,numel(perms))));
perms = randperm(numel(totalCreLM));
totalCreLM = totalCreLM(perms(1:min(n,numel(perms))));

save(fullfile(outputDirectory, 'RAW_DATA_aggregateActinSpeedVSdistanceToEdge_12um.mat'),...
  'totalConLP', 'totalConLM', 'totalCreLP', 'totalCreLM');

Y = [muTotalConLP, muTotalCreLP; muTotalConLM, muTotalCreLM];
S = [stdTotalConLP, stdTotalCreLP; stdTotalConLM, stdTotalCreLM];
N = [nTotalConLP, nTotalCreLP; nTotalConLM, nTotalCreLM];
SE = S ./ sqrt(N); % Standard Error of the mean

% Display bars
hFig = figure('Visible', 'off');
h = bar(Y, 'grouped');
set(gca, 'FontName', 'Helvetica', 'FontSize', 14);
set(gca, 'XTickLabel', {'0-1700 nm', '> 1700 nm'});
ylabel('Actin Speed Outside Adhesions (nm/min)');
hold on;

hAxes = get(h, 'Children');
XData1 = get(hAxes{1}, 'XData');
XData1 = mean(XData1(1:2:end, :), 1);
XData2 = get(hAxes{2}, 'XData');
XData2 = mean(XData2(1:2:end, :), 1);
X = [XData1', XData2'];

% Display Standard Error of the Mean
errorbar(reshape(X,4,1), reshape(Y,4,1), reshape(1.96 * SE,4,1),'xk');

% Perform the T-test on 2 control bins
pValue = ranksum(totalConLP, totalConLM);

if pValue < alpha
    y = max(Y(:));
    line(X(:,1), repmat(y + 80, 1, 2) + 10, 'Color', 'k', 'LineWidth', 2);
    text(mean(X(:,1)), y + 95, '*', 'FontName', 'Helvetica', 'FontSize', 24);    
end

% Perform the T-test on 2 null bins
pValue = ranksum(totalCreLP, totalCreLM);

if pValue < alpha
    y = max(Y(:));
    line(X(:,2), repmat(y + 130, 1, 2) + 10, 'Color', 'k', 'LineWidth', 2);
    text(mean(X(:,2)), y + 145, '*', 'FontName', 'Helvetica', 'FontSize', 24);    
end

% Perform the T-test on the 1st control and null bins
pValue = ranksum(totalConLP, totalCreLP);

if pValue < alpha
    y = max(Y(:));
    line(X(1,:), repmat(y + 20, 1, 2) + 10, 'Color', 'k', 'LineWidth', 2);
    text(mean(X(1,:)), y + 35, '*', 'FontName', 'Helvetica', 'FontSize', 24);    
end

% Perform the T-test on the 2nd control and null bins
pValue = ranksum(totalConLM, totalCreLM);
if pValue < alpha
    y = max(Y(:));
    line(X(2,:), repmat(y + 20, 1, 2) + 10, 'Color', 'k', 'LineWidth', 2);
    text(mean(X(2,:)), y + 35, '*', 'FontName', 'Helvetica', 'FontSize', 24);    
end

% Saving
fileName = fullfile(outputDirectory, 'aggregateActinSpeedOutsideAdhVSdistanceToCellEdge_fig4C.eps');
print(hFig, '-depsc', fileName);
fixEpsFile(fileName);
close(hFig);
