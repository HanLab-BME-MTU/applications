
%% Common variables between all plots

%define condition ID ('MM','MP','PM' or 'MM' - first M/P refers to VEGF, second MP refers to AAL)
condID = 'PP';

%get proper resSummaryInd and dataset names
eval(['resSummaryIndUse = resSummaryInd' condID ';'])
eval(['condName = condName' condID ';'])

%number of individual datasets
numDS = length(resSummaryIndUse);

%time list
timeList = cell(numDS,1);
for iDS = 1 : numDS
    timeList{iDS} = resSummaryIndUse(iDS).timeList(:,1);
end

%dataset colors
condColor = {[0 0 0],[0 0 1],[1 0 0],[0 1 0],[1 0 1],[0 1 1],[1 1 0],[0.7 0.5 0]};

% %try to match colors and dates - ANOTHER TIME
% %dates: 141028, 141030, 141031, 150128, 150129, 150224/25, 150302/03,
% %150307/10, 150311, 150312, 150317/18, 150319
% condColor = {[0 0 0],[0 0 1],[1 0 0],[0 1 0],[1 0 1],[0 1 1],[1 1 0],[0.7 0.5 0],[0.7 0.7 0.7]};

%flag whether to shift negative time so that first time point for all
%datasets is 0
shiftNegTime = 1;

%directory for saving figures
dir2save = '/project/biophysics/jaqaman_lab/vegf_tsp1/slee/VEGFR2/analysisCombinedDatasets/150320_analysis/figuresInd';

%% Absolute number of molecules in various motion classes

tmp = cell(numDS,1);
for iDS = 1 : numDS
    tmp{iDS} = resSummaryIndUse(iDS).numAbsClass;
end

figNameList = {'numAbsImm','numAbsConf','numAbsFree','numAbsDir',...
    'numAbsUndet','numAbsDet','numAbsTot'};
for iFig = 1 : length(figNameList)
    figNameList{iFig} = [figNameList{iFig} condID];
end

yaxisUnits = '(molecules)';

plotResMultipleDatasets(tmp,timeList,condName,condColor,figNameList,dir2save,yaxisUnits,shiftNegTime);

%% Normalized number of molecules in various motion classes

tmp = cell(numDS,1);
for iDS = 1 : numDS
    tmp{iDS} = resSummaryIndUse(iDS).numNorm0Class;
end

figNameList = {'numNorm0Imm','numNorm0Conf','numNorm0Free','numNorm0Dir',...
    'numNorm0Undet','numNorm0Det','numNorm0Tot'};
for iFig = 1 : length(figNameList)
    figNameList{iFig} = [figNameList{iFig} condID];
end

yaxisUnits = '(unitless)';

plotResMultipleDatasets(tmp,timeList,condName,condColor,figNameList,dir2save,yaxisUnits,shiftNegTime);

%% Probability of various motion classes

tmp = cell(numDS,1);
for iDS = 1 : numDS
    tmp{iDS} = resSummaryIndUse(iDS).probClass;
end

figNameList = {'probImm','probConf','probFree','probDir','probDet'};
for iFig = 1 : length(figNameList)
    figNameList{iFig} = [figNameList{iFig} condID];
end

yaxisUnits = '(unitless)';

plotResMultipleDatasets(tmp,timeList,condName,condColor,figNameList,dir2save,yaxisUnits,shiftNegTime);

%% Diffusion coefficient in various motion classes

tmp = cell(numDS,1);
for iDS = 1 : numDS
    tmp{iDS} = resSummaryIndUse(iDS).diffCoefClass;
end

figNameList = {'diffCoefImm','diffCoefConf','diffCoefFree','diffCoefDir','diffCoefUndet'};
for iFig = 1 : length(figNameList)
    figNameList{iFig} = [figNameList{iFig} condID];
end

yaxisUnits = '(pixels^2/frame)';

plotResMultipleDatasets(tmp,timeList,condName,condColor,figNameList,dir2save,yaxisUnits,shiftNegTime);

%% Confinement radius in various motion classes

tmp = cell(numDS,1);
for iDS = 1 : numDS
    tmp{iDS} = resSummaryIndUse(iDS).confRadClass;
end

figNameList = {'confRadImm','confRadConf','confRadFree','confRadDir','confRadUndet'};
for iFig = 1 : length(figNameList)
    figNameList{iFig} = [figNameList{iFig} condID];
end

yaxisUnits = '(pixels)';

plotResMultipleDatasets(tmp,timeList,condName,condColor,figNameList,dir2save,yaxisUnits,shiftNegTime);

%% amplitude in various motion classes

tmp = cell(numDS,1);
for iDS = 1 : numDS
    tmp{iDS} = resSummaryIndUse(iDS).ampClass;
end

figNameList = {'ampImm','ampConf','ampFree','ampDir','ampUndet'};
for iFig = 1 : length(figNameList)
    figNameList{iFig} = [figNameList{iFig} condID];
end

yaxisUnits = '(arbitrary units)';

plotResMultipleDatasets(tmp,timeList,condName,condColor,figNameList,dir2save,yaxisUnits,shiftNegTime);

%% amplitude in first and last 20 frames

tmp = cell(numDS,1);
for iDS = 1 : numDS
    tmp{iDS} = resSummaryIndUse(iDS).ampFL20;
end

figNameList = {'ampF20','ampL20'};
for iFig = 1 : length(figNameList)
    figNameList{iFig} = [figNameList{iFig} condID];
end

yaxisUnits = '(arbitrary units)';

plotResMultipleDatasets(tmp,timeList,condName,condColor,figNameList,dir2save,yaxisUnits,shiftNegTime);

%% rates of merging and splitting

tmp = cell(numDS,1);
for iDS = 1 : numDS
    tmp{iDS} = resSummaryIndUse(iDS).rateMS;
end

figNameList = {'rateMerge','rateSplit'};
for iFig = 1 : length(figNameList)
    figNameList{iFig} = [figNameList{iFig} condID];
end

yaxisUnits = '(per frame)';

plotResMultipleDatasets(tmp,timeList,condName,condColor,figNameList,dir2save,yaxisUnits,shiftNegTime);

