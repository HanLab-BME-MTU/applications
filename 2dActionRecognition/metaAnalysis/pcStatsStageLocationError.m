%% Statistics on stage location errors across all experiments
function [] = pcStatsStageLocationError()

addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition/metaAnalysis/'));

close all;

analysisDirname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/';
metaDataFname = [analysisDirname 'MetaData/ExperimentsAll20170509.mat'];

load(metaDataFname);

translations = [];
for i = 1 : metaData.tasks.N
    curExp = metaData.tasks.exps(i);
    curTask = metaData.tasks.tasks(i);
    curFname = metaData.experiments.fnames{curExp};    
    
    fname = [analysisDirname 'All/' curFname filesep...
        curFname '_s' sprintf('%02d',curTask) filesep...
        'MF' filesep 'mf' filesep 'stageLocationError.mat'];
    
    if ~exist(fname,'file')
        fname = [analysisDirname 'All/' curFname filesep...
            curFname '_s' sprintf('%d',curTask) filesep...
            'MF' filesep 'mf' filesep 'stageLocationError.mat'];
    end
    
    if ~exist(fname,'file')
        continue;
    end
    
    load(fname);
    translations = [translations stageShift];    
end

bins = 0.25:0.5:3;
[nelements, ~] = hist(translations,bins);
transDist = nelements ./ sum(nelements);
h = figure;
hold on;
bar(bins,transDist,'r');
xlabel(xlabelStr,'FontSize',22);
ylabel('Fraction','FontSize',22);
haxes = get(h,'CurrentAxes');
set(haxes,'XLim',xlimVals);
set(haxes,'XTick',xlimVals(1):xlimStep:xlimVals(2));
set(haxes,'XTickLabel',xlimVals(1):xlimStep:xlimVals(2));
set(haxes,'YLim',ylimVals);
set(haxes,'YTick',ylimVals(1):ylimStep:ylimVals(2));
set(haxes,'YTickLabel',ylimVals(1):ylimStep:ylimVals(2));
set(haxes,'FontSize',22);
set(h,'Color','none');
hold off;
export_fig(outFname);

end