function [groupedData] = groupArmaCoefficients(data,options)
%GROUPARMACOEFFICIENTS groups models according to similarity of their ARMA coefficients
%
% SYNOPSIS: [out] = groupArmaCoefficients(data,options)
%
% INPUT
%
% OUTPUT
%
% REMARKS
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 22-May-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%============================
%% Test Input / Load Data
%============================

groupedData = [];
%def_mode = 4;

% default cutoffs. Negative integers are number of groups, postitive reals
% are probability cutoffs. Now that -log10(prob) is being used for
% clustering, you have to specify the log of the cutoff probability
def_options.wnv1 = -2;
def_options.arma = -2;%-log10(5e-5); %-2
def_options.wnv2 = -log10(1e-12);

% tbd: test input, check options, set defaults

if nargin == 0 || isempty(data)
    % load data. have user select list of strains
    [strainFile, dataPath] = uigetfile('strains_*.mat','Load strain list!');
    if any(strainFile == 0)
        disp('--groupArmaCoefficients aborted')
        return
    end
    
    % read strainInfo so that user may select the files to be analyzed
    % directly without having to wait for the data to be read from the
    % server
    strainInfo = load(fullfile(dataPath,strainFile));
    
     % strainInfo is a structure in itself
    fn = fieldnames(strainInfo);
    strainInfo = strainInfo.(fn{1});


    % let the user select the strains
    [nameList{1:length(strainInfo)}] = deal(strainInfo.name);
    selectionIdx = listSelectGUI(nameList);

    if isempty(selectionIdx)
        disp('--groupArmaCoefficients aborted')
        return
    end
    
    
   
    
    % read armaData, lengthData. StrainInfo is a structure with
    % the dataFile-names, while armaData and lengthData are collections of
    % files that are labeled according to the dataFile-names.
     dataSuffix = regexp(strainFile,'strains_([\w]+).mat','tokens');
    dataSuffix = char(dataSuffix{1});
    % armaData is a collection of many files, while strainInfo is already
    % the correct structure that contains the optimal model for each
    armaData = load(fullfile(dataPath,sprintf('resFitVelAndLen_%s',dataSuffix)));
    
    % try load length file
    try
        lengthData = load(fullfile(dataPath,sprintf('lengthSeries_%s',dataSuffix)));
    catch
        lengthData = [];
    end

   


    % get number of data
    nData = length(selectionIdx);


    % read the correct model. If orderLen is [], skip.
    for i=nData:-1:1,
        if isempty(strainInfo(selectionIdx(i)).orderLen)
            % remove only if there is data already
            if exist('data','var')
                data(i) = [];
                goodIdx(i) = [];
            end
        else
            % reorder data fields, in case there is a change in field-order in the future
            %             if exist('data','var')
            %             data = orderfields(armaData.(sprintf('fitLen%s',strainInfo(selectionIdx(i)).name))...
            %             (strainInfo(selectionIdx(i)).orderLen(1)+1,strainInfo(selectionIdx(i)).orderLen(2)+1),data);
            %             end
            data(i)=...
                armaData.(sprintf('fitLen%s',strainInfo(selectionIdx(i)).name))...
                (strainInfo(selectionIdx(i)).orderLen(1)+1,strainInfo(selectionIdx(i)).orderLen(2)+1);
            goodIdx(i) = true;
        end
    end
    % reset nData
    data = data(:);
    nData = length(data);


    % read strainInfo into data
    [data.name] = deal(strainInfo(selectionIdx(goodIdx)).name);
    [data.orderLen] = deal(strainInfo(selectionIdx(goodIdx)).orderLen);
    
    % read length into data (if there is any). We can't do this above,
    % because filling up the structure requires the fields to be the same
    if ~isempty(lengthData)
        % we already know what part of the selection is good
        for i=selectionIdx(goodIdx)'
            data(i).lengthSeries = ...
                lengthData.(sprintf('length%s',strainInfo(selectionIdx(i)).name));
        end
    end

end

if nargin < 2 || isempty(options)
    options = struct;
end
% assign default options
defaultOptions = fieldnames(def_options);
for fn = 1:length(defaultOptions)
    if ~isfield(options,defaultOptions{fn}) || isempty(options.(defaultOptions{fn}))
        options.(defaultOptions{fn}) = def_options.(defaultOptions{fn});
    end
end

%============================


%============================
%% Group according to WNV
%============================

groupingOptions.cutoff = options.wnv1;
groupingOptions.labels = strvcat(data.name);
[links, groupIdx, groupedData, linkData] = groupData(data,'groupArma_distance_WNV',groupingOptions);

% label axis
set(get(groupedData.plotHandles.axesH,'XLabel'),'String','-^1^0log(probability)')
set(groupedData.plotHandles.figureH,'Name','WNV - Round 1');


%===========================


%===========================
%% Group according to ARMA
%===========================

%options.mode = 2;
groupingOptions.cutoff = options.arma;

nGroups = length(groupedData.collectedData);

% plot
plotRows = ceil(sqrt(nGroups));
plotCols = ceil(nGroups/plotRows);
fh = figure('Name','ARMA-clustering');


for iGroup = 1:nGroups
    figure(fh);
    sh = subplot(plotRows,plotCols,iGroup);
    % group, find subgroups
    groupingOptions.labels = strvcat(groupedData.collectedData(iGroup).data.name);
    [groupedData.collectedData(iGroup).links,groupedData.collectedData(iGroup).groupIdx,...
        groupedData.collectedData(iGroup).subGroups,groupedData.collectedData(iGroup).linkData] = groupData(...
        groupedData.collectedData(iGroup).data,'groupArma_distance_ARMA',groupingOptions);

    % steal all we can from the figure
    ah = groupedData.collectedData(iGroup).subGroups.plotHandles.axesH;
    % steal lines
    copyobj(get(ah,'Children'),sh)
    % steal y-labels
    set(sh,'YTick',get(ah,'YTick'),'YTickLabel',get(ah,'YTickLabel'),'YLim',get(ah,'YLim'))
    set(get(sh,'XLabel'),'String','-^1^0log(probability)')

    % close figure
    delete(groupedData.collectedData(iGroup).subGroups.plotHandles.figureH);

end

%===========================
%% Group according to WNV
%===========================

%options.mode = 2;
groupingOptions.cutoff = options.wnv2;


for iGroup = 1:nGroups

    nSubGroups = length(groupedData.collectedData(iGroup).subGroups.collectedData);

    % plot
    plotRows = ceil(sqrt(nSubGroups));
    plotCols = ceil(nSubGroups/plotRows);
    fh = figure('Name',sprintf('WNV - Round 2. Arma-Group %i',iGroup));


    for jGroup = 1:nSubGroups
        figure(fh);
        sh = subplot(plotRows,plotCols,jGroup);
        % group
        if length(groupedData.collectedData(iGroup).subGroups.collectedData(jGroup).data)>1
            groupingOptions.labels = strvcat(groupedData.collectedData(iGroup).subGroups.collectedData(jGroup).data.name);
            [groupedData.collectedData(iGroup).subGroups.collectedData(jGroup).links,...
                groupedData.collectedData(iGroup).subGroups.collectedData(jGroup).groupIdx,...
                groupedData.collectedData(iGroup).subGroups.collectedData(jGroup).subGroups,...
                groupedData.collectedData(iGroup).subGroups.collectedData(jGroup).linkData] = groupData(...
                groupedData.collectedData(iGroup).subGroups.collectedData(jGroup).data,'groupArma_distance_WNV',groupingOptions);
            xlabel('-^1^0log(probability)')
        end

        if isfield(groupedData.collectedData(iGroup).subGroups.collectedData(jGroup),'subGroups') &&...
                ~isempty(groupedData.collectedData(iGroup).subGroups.collectedData(jGroup).subGroups)
            % steal all we can from the figure
            ah = groupedData.collectedData(iGroup).subGroups.collectedData(jGroup).subGroups.plotHandles.axesH;
            % steal lines
            copyobj(get(ah,'Children'),sh)
            % steal y-labels
            set(sh,'YTick',get(ah,'YTick'),'YTickLabel',get(ah,'YTickLabel'),'YLim',get(ah,'YLim'))
            set(get(sh,'XLabel'),'String','-^1^0log(probability)')

            % close figure
            delete(groupedData.collectedData(iGroup).subGroups.collectedData(jGroup).subGroups.plotHandles.figureH);
        end
    end
end