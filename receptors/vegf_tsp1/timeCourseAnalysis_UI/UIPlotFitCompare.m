function [] = UIPlotFitCompare()
%plots the point by point p-value comparrison between two fitted curves from timeCourseAnalysis 
%
%Tae H Kim, July 2015

%% User Input
%get figure data
[fileName, filePath] = uigetfile('*.mat', 'Select figureData', 'MultiSelect', 'off');
outputDir = [filePath 'figure_p'];
%loads figure data
load([filePath fileName]);
%ask which figures' p-value comparrisons should be plotted
stringList = arrayfun(@(x) [x.titleBase ' ' x.titleVariable], figureData, 'UniformOutput', false);
figIndx = listdlg('PromptString','Select Movie:', 'SelectionMode','multiple', 'ListString', stringList, 'ListSize', [300, 400]);
%make directory where new plots will be saved
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% Initialize
nCond = numel(commonInfo.times);
nFig = numel(figIndx);
%obtain condition pairs
nPair = nCond * (nCond - 1) / 2;
pairIndx = cell(1, nPair);
pairName = cell(1, nPair);
iPair = 0;
for indx1 = 1:nCond-1
    for indx2 = indx1+1:nCond
        iPair = iPair + 1;
        pairIndx{iPair} = [indx1, indx2];
        pairName{iPair} = [commonInfo.conditions{indx1} ' vs ' commonInfo.conditions{indx2}];
    end
end
%colors
condColorAll = {...
    [0 0 0],... %black 1
    [0 0 1],... %blue 2
    [0 0.5 0],... %dark green 3
    [1 0 0],... %red 4
    [1 0.6 0.3],... %orange 5
    [0 1 0],... %green 6
    [1 0 1],... %magenta 7
    [0 1 1],... %cyan 8
    [1 1 0],... %yellow 9
    [0.7 0.5 0],... %brown 10
    [0.7 0.7 0.7]... %gray 11
    [0.5 0.5 1],... %purple 12
    [0.3 0.8 1],... %light blue 13
    };
colors = condColorAll(mod(1:nPair, 13) + 1);

%% Plot
%Progress Display
progressText_Increment('Ploting figures', nFig);
%Plots selected figures
lineObj(nPair) = line();
for iFig = figIndx
    hold off
    figObj = figure;
    hold on;
    for iPair = 1:nPair
        indx1 = pairIndx{iPair}(1);
        indx2 = pairIndx{iPair}(2);
        timeIndx = figureData(iFig).fitCompare(indx1, indx2).timeIndx;
        lineObj(iPair) = plot(commonInfo.compareTimes{timeIndx}', figureData(iFig).fitCompare(indx1, indx2).p', 'color', colors{iPair});
    end
    legend(lineObj, pairName);
    title([figureData(iFig).titleBase, ' ', figureData(iFig).titleVariable, ' p-values between curves']);
    xlabel('Time (min)');
    ylabel('p');
    xlim manual;
    ylim([0, 0.05]);
    savefig(figObj, [outputDir, filesep, figureData(iFig).titleBase, ' ', figureData(iFig).titleVariable, '.fig']);
    close(figObj);
    %progress display
    progressText_Increment();
end
end

