function aggregateActinSpeedVSadhesionLength(conDirectory, creDirectory, outputDirectory, alpha, n)

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
nConNasc = zeros(nConMovies,1);
nConElong = zeros(nConMovies,1);
nCreNasc = zeros(nCreMovies,1);
nCreElong = zeros(nCreMovies,1);

% Average Actin speed in bin 1 and 2 for control and null movies
muConNasc = zeros(nConMovies,1);
muConElong = zeros(nConMovies,1);
muCreNasc = zeros(nCreMovies,1);
muCreElong = zeros(nCreMovies,1);

% Standard deviation of Actin speed in bin 1 and 2 for control and null
% movies
stdConNasc = zeros(nConMovies,1);
stdConElong = zeros(nConMovies,1);
stdCreNasc = zeros(nCreMovies,1);
stdCreElong = zeros(nCreMovies,1);

totalConNasc = cell(nConNasc,1);
totalConElong = cell(nConElong,1);
totalCreNasc = cell(nCreNasc,1);
totalCreElong = cell(nCreElong,1);

% Nasc in controls
for iMovie = 1:nConMovies    
    % load data
    load(conPaths{iMovie});
    
    totalConNasc{iMovie} = data{1};
    totalConElong{iMovie} = data{2};

    nConNasc(iMovie) = numel(data{1}); %#ok<*USENS>
    muConNasc(iMovie)= mean(data{1});
    stdConNasc(iMovie) = std(data{1});
    
    nConElong(iMovie) = numel(data{2});
    muConElong(iMovie) = mean(data{2});
    stdConElong(iMovie) = std(data{2});
    
    fprintf(1, '%s: (BIN 1) number of adhesion counts = %d, mean = %f, std = %f\n', getDirFromPath(conPaths{iMovie}), nConNasc(iMovie), muConNasc(iMovie), stdConNasc(iMovie));
    fprintf(1, '%s: (BIN 2) number of adhesion counts = %d, mean = %f, std = %f\n', getDirFromPath(conPaths{iMovie}), nConElong(iMovie), muConElong(iMovie), stdConElong(iMovie));
end

% Nasc in nulls
for iMovie = 1:nCreMovies    
    % load data
    load(crePaths{iMovie});

    totalCreNasc{iMovie} = data{1};
    totalCreElong{iMovie} = data{2};

    nCreNasc(iMovie) = numel(data{1}); %#ok<*USENS>
    muCreNasc(iMovie)= mean(data{1});
    stdCreNasc(iMovie) = std(data{1});
    
    nCreElong(iMovie) = numel(data{2});
    muCreElong(iMovie) = mean(data{2});
    stdCreElong(iMovie) = std(data{2});

    fprintf(1, '%s: (BIN 1) number of adhesion counts = %d, mean = %f, std = %f\n', getDirFromPath(crePaths{iMovie}), nCreNasc(iMovie), muCreNasc(iMovie), stdCreNasc(iMovie));
    fprintf(1, '%s: (BIN 2) number of adhesion counts = %d, mean = %f, std = %f\n', getDirFromPath(crePaths{iMovie}), nCreElong(iMovie), muCreElong(iMovie), stdCreElong(iMovie));
end

nTotalConNasc = sum(nConNasc);
nTotalConElong = sum(nConElong);
nTotalCreNasc = sum(nCreNasc);
nTotalCreElong = sum(nCreElong);

muTotalConNasc = mean(muConNasc);
muTotalConElong = mean(muConElong);
stdTotalConNasc = sqrt((1/nConMovies) * sum(stdConNasc.^2));
stdTotalConElong = sqrt((1/nConMovies) * sum(stdConElong.^2));

muTotalCreNasc = mean(muCreNasc);
muTotalCreElong = mean(muCreElong);
stdTotalCreNasc = sqrt((1/nCreMovies) * sum(stdCreNasc.^2));
stdTotalCreElong = sqrt((1/nCreMovies) * sum(stdCreElong.^2));

totalConNasc = vertcat(totalConNasc{:});
totalConElong = vertcat(totalConElong{:});
totalCreNasc = vertcat(totalCreNasc{:});
totalCreElong = vertcat(totalCreElong{:});

% sub-sample data:
perms = randperm(numel(totalConNasc));
totalConNasc = totalConNasc(perms(1:min(n,numel(perms))));
perms = randperm(numel(totalConElong));
totalConElong = totalConElong(perms(1:min(n,numel(perms))));
perms = randperm(numel(totalCreNasc));
totalCreNasc = totalCreNasc(perms(1:min(n,numel(perms))));
perms = randperm(numel(totalCreElong));
totalCreElong = totalCreElong(perms(1:min(n,numel(perms))));

save(fullfile(outputDirectory, 'RAW_DATA_aggregateActinSpeedVSadhesionLength_12um.mat'),...
  'totalConNasc', 'totalConElong', 'totalCreNasc', 'totalCreElong');

Y = [muTotalConNasc, muTotalCreNasc; muTotalConElong, muTotalCreElong];
S = [stdTotalConNasc, stdTotalCreNasc; stdTotalConElong, stdTotalCreElong];
N = [nTotalConNasc, nTotalCreNasc; nTotalConElong, nTotalCreElong];
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

% % Perform the U-test on 2 control bins
A = randperm(numel(totalConNasc));
B = randperm(numel(totalConElong));
pValue = ranksum(totalConNasc,totalConElong);

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

% Perform the U-test on 2 null bins
pValue = ranksum(totalCreNasc,totalCreElong);

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

% Perform the U-test on the 1st control and null bins
pValue = ranksum(totalConNasc, totalCreNasc);

if pValue < alpha
    y = max(Y(:));
    line(X(1,:), repmat(y + 20, 1, 2) + 10, 'Color', 'k', 'LineWidth', 2);
    text(mean(X(1,:)), y + 35, '*', 'FontName', 'Helvetica', 'FontSize', 24);    
end

% Perform the U-test on the 2nd control and null bins
pValue = ranksum(totalConElong, totalCreElong);

if pValue < alpha
    y = max(Y(:));
    line(X(2,:), repmat(y + 20, 1, 2) + 10, 'Color', 'k', 'LineWidth', 2);
    text(mean(X(2,:)), y + 35, '*', 'FontName', 'Helvetica', 'FontSize', 24);    
end

% Saving
fileName = fullfile(outputDirectory, 'aggregateActinSpeedVSadhesionLength_fig4C_Nascon2um_Elongon12um.eps');
print(hFig, '-depsc', fileName);
fixEpsFile(fileName);
close(hFig);
