function aggregateActinSpeedVSadhesionLength(conDirectory, creDirectory, outputDirectory, alpha)

cmd1 = ['find ' conDirectory ' -name ''*_DATA_12000nm.mat'''];
cmd2 = ['find ' creDirectory ' -name ''*_DATA_12000nm.mat'''];

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

% Number of adhesion counts in bin 1 and 2 for control and null movies
nConBin1 = zeros(nConMovies,1);
nConBin2 = zeros(nConMovies,1);
nCreBin1 = zeros(nCreMovies,1);
nCreBin2 = zeros(nCreMovies,1);

% Average Actin speed in bin 1 and 2 for control and null movies
muConBin1 = zeros(nConMovies,1);
muConBin2 = zeros(nConMovies,1);
muCreBin1 = zeros(nCreMovies,1);
muCreBin2 = zeros(nCreMovies,1);

% Standard deviation of Actin speed in bin 1 and 2 for control and null
% movies
stdConBin1 = zeros(nConMovies,1);
stdConBin2 = zeros(nConMovies,1);
stdCreBin1 = zeros(nCreMovies,1);
stdCreBin2 = zeros(nCreMovies,1);

% Bin1 in controls
for iMovie = 1:nConMovies    
    % load data
    load(conPaths{iMovie});
    
    nConBin1(iMovie) = numel(data{1}); %#ok<*USENS>
    muConBin1(iMovie)= mean(data{1});
    stdConBin1(iMovie) = std(data{1});
    
    nConBin2(iMovie) = numel(data{2});
    muConBin2(iMovie) = mean(data{2});
    stdConBin2(iMovie) = std(data{2});
    
    fprintf(1, '%s: (BIN 1) number of adhesion counts = %d, mean = %f, std = %f\n', getDirFromPath(conPaths{iMovie}), nConBin1(iMovie), muConBin1(iMovie), stdConBin1(iMovie));
    fprintf(1, '%s: (BIN 2) number of adhesion counts = %d, mean = %f, std = %f\n', getDirFromPath(conPaths{iMovie}), nConBin2(iMovie), muConBin2(iMovie), stdConBin2(iMovie));
end

% Bin1 in nulls
for iMovie = 1:nCreMovies    
    % load data
    load(crePaths{iMovie});
    
    nCreBin1(iMovie) = numel(data{1}); %#ok<*USENS>
    muCreBin1(iMovie)= mean(data{1});
    stdCreBin1(iMovie) = std(data{1});
    
    nCreBin2(iMovie) = numel(data{2});
    muCreBin2(iMovie) = mean(data{2});
    stdCreBin2(iMovie) = std(data{2});

    fprintf(1, '%s: (BIN 1) number of adhesion counts = %d, mean = %f, std = %f\n', getDirFromPath(crePaths{iMovie}), nCreBin1(iMovie), muCreBin1(iMovie), stdCreBin1(iMovie));
    fprintf(1, '%s: (BIN 2) number of adhesion counts = %d, mean = %f, std = %f\n', getDirFromPath(crePaths{iMovie}), nCreBin2(iMovie), muCreBin2(iMovie), stdCreBin2(iMovie));
end

nTotalConBin1 = sum(nConBin1);
nTotalConBin2 = sum(nConBin2);
nTotalCreBin1 = sum(nCreBin1);
nTotalCreBin2 = sum(nCreBin2);

muTotalConBin1 = mean(muConBin1);
muTotalConBin2 = mean(muConBin2);
stdTotalConBin1 = sqrt((1/nConMovies) * sum(stdConBin1.^2));
stdTotalConBin2 = sqrt((1/nConMovies) * sum(stdConBin2.^2));

muTotalCreBin1 = mean(muCreBin1);
muTotalCreBin2 = mean(muCreBin2);
stdTotalCreBin1 = sqrt((1/nCreMovies) * sum(stdCreBin1.^2));
stdTotalCreBin2 = sqrt((1/nCreMovies) * sum(stdCreBin2.^2));

Y = [muTotalConBin1, muTotalCreBin1; muTotalConBin2, muTotalCreBin2];
S = [stdTotalConBin1, stdTotalCreBin1; stdTotalConBin2, stdTotalCreBin2];
N = [nTotalConBin1, nTotalCreBin1; nTotalConBin2, nTotalCreBin2];
SE = S ./ sqrt(N); % Standard Error of the mean

% Display bars
hFig = figure('Visible', 'off');
h = bar(Y, 'grouped');
set(gca, 'FontName', 'Helvetica', 'FontSize', 14);
set(gca, 'XTickLabel', {'Diffraction Limited', 'Elongating'});
ylabel('Actin Speed (nm/min)');
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
n1 = nTotalConBin1;
n2 = nTotalConBin2;
m1 = muTotalConBin1;
m2 = muTotalConBin2;
v1 = stdTotalConBin1^2;
v2 = stdTotalConBin2^2;
T = (m1 - m2) / sqrt((((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2)) * (1/n1 + 1/n2));
pValue = (1 - tcdf(T, n1 + n2 - 2));

if pValue < alpha
    y = max(Y(:));
    line(X(:,1), repmat(y + 80, 1, 2) + 10, 'Color', 'k', 'LineWidth', 2);
    if pValue
        str = sprintf('* (p-value=%E)', pValue);
    else
        str = sprintf('* (p-value=0)');
    end
    text(mean(X(:,1)), y + 110, str, 'FontName', 'Helvetica', 'FontSize', 12);    
end

% Perform the T-test on 2 null bins
n1 = nTotalCreBin1;
n2 = nTotalCreBin2;
m1 = muTotalCreBin1;
m2 = muTotalCreBin2;
v1 = stdTotalCreBin1^2;
v2 = stdTotalCreBin2^2;
T = (m1 - m2) / sqrt((((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2)) * (1/n1 + 1/n2));
pValue = (1 - tcdf(T, n1 + n2 - 2));

if pValue < alpha
    y = max(Y(:));
    line(X(:,2), repmat(y + 130, 1, 2) + 10, 'Color', 'k', 'LineWidth', 2);
    if pValue
        str = sprintf('* (p-value=%E)', pValue);
    else
        str = sprintf('* (p-value=0)');
    end
    text(mean(X(:,1)), y + 160, str, 'FontName', 'Helvetica', 'FontSize', 12);    
end

% Perform the T-test on the 1st control and null bins
n1 = nTotalCreBin1;
n2 = nTotalConBin1;
m1 = muTotalCreBin1;
m2 = muTotalConBin1;
v1 = stdTotalCreBin1;
v2 = stdTotalConBin1;

T = (m1 - m2) / sqrt((((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2)) * (1/n1 + 1/n2));
pValue = (1 - tcdf(T, n1 + n2 - 2));

if pValue < alpha
    y = max(Y(:));
    line(X(1,:), repmat(y + 20, 1, 2) + 10, 'Color', 'k', 'LineWidth', 2);
    text(mean(X(1,:)), y + 35, '*', 'FontName', 'Helvetica', 'FontSize', 24);    
end

% Perform the T-test on the 2nd control and null bins
n1 = nTotalCreBin2;
n2 = nTotalConBin2;
m1 = muTotalCreBin2;
m2 = muTotalConBin2;
v1 = stdTotalCreBin2;
v2 = stdTotalConBin2;

T = (m1 - m2) / sqrt((((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2)) * (1/n1 + 1/n2));
pValue = (1 - tcdf(T, n1 + n2 - 2));

if pValue < alpha
    y = max(Y(:));
    line(X(2,:), repmat(y + 20, 1, 2) + 10, 'Color', 'k', 'LineWidth', 2);
    text(mean(X(2,:)), y + 35, '*', 'FontName', 'Helvetica', 'FontSize', 24);    
end

% Saving
fileName = fullfile(outputDirectory, 'aggregateActinSpeedVSadhesionLength_fig4C_bin1on2um_bin2on12um.eps');
print(hFig, '-depsc', fileName);
fixEpsFile(fileName);
close(hFig);
