function [groupIndex,groupedData] = groupArmaDescriptors(data,options)
%GROUPARMADESCRIPTORS groups models according to similarity of their ARMA coefficients
%
% SYNOPSIS: groupedData = groupArmaDescriptors(data,options)
%
% INPUT     data (opt): structure similar to the output of armaxFitKalman, with
%                 two additional fields: orderLen (number of fit
%                 parameters), name (name of the dataSet). If empty, the
%                 code will ask for a file called strain_???.mat, an then
%                 look for a file called resFitVelAndLen_???.mat, and,
%                 possibly, lengthSeries_???.mat
%           options (opt): Structure with options for grouping data sets
%                 (### = wnv1, arma, or wnv2 for the first round of
%                 clustering according to white noise variance, clustering
%                 according to arma coefficients, and the second round of
%                 clustering according to white noise variance,
%                 respectively.
%                 ###_cutoff - if negative integer: fixed number of groups
%                              into which to split data.
%                              if positive real number (0...1): p-value
%                              indicating significant differences
%                 ###_mode -   two-element vector describing how individual
%                              data sets should be combined
%                              [0,0] - no recalculation of ARMA descriptors
%                              [1,p] - recalculation until a threshold p is reached
%                              [2,0] - recalculation for all groups
%                              [3,p,a] - recalculation until a threshold p
%                              is reached (p==0 => recalc for all groups).
%                              Use only a subset of movies from each group,
%                              so that the total number of observations of
%                              the recalculated group is the average of the
%                              number of the observations of all the
%                              combined sets. a==1 chooses roughly the same
%                              number of timepoints from each set, a==2
%                              takes timepoints according to the size of
%                              the data sets.
%                 multiply     When grouping individual trajectories, set
%                              this option to a positive integer n. The
%                              individual trajectories will be repeated n
%                              times to allow recalculation with
%                              subsampling. Default: 0
%
% OUTPUT    groupIndex : nData-by-3 array that indicates for every data set
%                        how it was grouped at levels 1,2, and 3.
%           groupedData : output structure with fields
%                .collectedData(1:n) groups of data as divided according to
%                    options.###_cutoff with fields
%                   .data : original data sets belonging to this group
%                   .dataIdx : index of data into input data order
%                   .links : links between data sets within the group (as
%                            returned by linkage.m)
%                   .linkData : additional data describing the links. Cols
%                       1: separate groups according to cutoff?
%                       2: -log10(prob)
%                       3: Fratio
%                       4: n1
%                       5: n2
%                       6: minLogProb between individual sets
%                       7: maxLogProb between individual sets
%                       8: 1 if extrapolated log
%                   .groupIdx : nSets-by-(nSets-1) array describing the
%                       coalescence of data sets into groups
%                   .subGroups : groupedData of subgroups as determined by
%                       cutoff
%                .collectedGrouping : links array of the relation between
%                    the n groups of collectedData
%                .plotHandles : structure with handles to the figures
%                   .figureH : figure handle
%                   .axesH   : axes handle
%                   .lineH   : handles to the lines of the dendrogram
%
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
if nargout > 0
    groupedData = [];
end

% default cutoffs. Negative integers are number of groups, postitive reals
% are probability cutoffs.
% mode describes whether the individual data sets should be combined, and
% to what level of difference.
% [0,0] - no recalculation of ARMA descriptors
% [1,p] - recalculation until a threshold p is reached
% [2,0] - recalculation for all groups
% [3,p,a] - recalculation until a threshold p is reached (p==0 => recalc
% for all groups). Use only a subset of movies from each group, so that the
% total number of observations of the recalculated group is the average of
% the number of the observations of all the combined sets. a==1 chooses
% roughly the same number of timepoints from each set, a==2 takes
% timepoints according to the size of the data sets.
def_options.wnv1_cutoff = -1;
def_options.wnv1_mode = [0 1e-12];
def_options.arma_cutoff = 5e-5;
def_options.arma_mode = [3,5e-5,1];%[0, 5e-5];
def_options.wnv2_cutoff = 1e-12;
def_options.wnv2_mode = [0 1e-12];
def_options.plot = 1; % plot results
def_options.multiply = 0; % don't repeat data
% remember system dialog setting
sysDialogState = getappdata(0,'UseNativeSystemDialogs');

% set options now, because we might need them for loading data
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

% tbd: test input, check options, set defaults

if nargin == 0 || isempty(data)
    % load data. have user select list of strains. Allow selection of
    % multiple files in the same directory (with while-loop like in
    % trajectoryAnalysis, it could be extended to span multiple
    % directories, or maybe Matlab imporves uigetfile)
    % 07/28/06 - there's a bug in Matlab/Linux. Since we're limited to one
    % directory at the moment, use the listSelectGUI to load strains
    % 08/10/06 - found workaround
    if isunix
        setappdata(0,'UseNativeSystemDialogs',false)
    end
    [strainFile, dataPath] = uigetfile('strains_*.mat','Load strain list(s)!','MultiSelect','on');
    if ~iscell(strainFile) && any(strainFile == 0)
        disp('--groupArmaCoefficients aborted')
        return
    end
    %     else
    %         dataPath = uigetdir('Coose directory of strain list(s)');
    %         strainFile = searchFiles('strains_','',dataPath,0);
    %         if isempty(strainFile)
    %             disp('--groupArmaCoefficients: no strains_* file found')
    %             return
    %         end
    %         % only remember the fileNames
    %         strainFile = strainFile(:,1);
    %
    %         % have the user choose
    %         [selectionIdx] = listSelectGUI(strainFile);
    %         if isempty(selectionIdx)
    %             disp('--groupArmaCoefficients aborted')
    %             return
    %         else
    %             strainFile = strainFile(selectionIdx);
    %         end
    %     end

    % make strainFile into a cell, so that we can use the code for
    % multiSelections with one selection.
    if ~iscell(strainFile)
        strainFile = {strainFile};
    end

    % read strainInfo so that user may select the files to be analyzed
    % directly without having to wait for the data to be read from the
    % server


    % load data in loop. This is, unfortunately, sensitive to changes in the
    % list and ordering of fields. One possible solution would be to
    % define the key list of fields above, and to only read those

    % load first strainInfo, then add the rest to the list, so that we
    % automatically get the proper list of fieldnames

    strainInfo = load(fullfile(dataPath,strainFile{1}));
    % strainInfo is a structure in itself
    fn = fieldnames(strainInfo);
    strainInfo = strainInfo.(fn{1});

    % loop for all the other files. If one file, the loop isn't executed
    for iFile = 2:length(strainFile)
        strainInfoTmp = load(fullfile(dataPath,strainFile{iFile}));
        fnTmp = fieldnames(strainInfoTmp);
        strainInfoTmp = strainInfoTmp.(fnTmp{1});
        % reorder fields, so that we're at least insensitive to that
        strainInfoTmp = orderfields(strainInfoTmp,strainInfo);
        strainInfo = [strainInfo(:);strainInfoTmp(:)];
    end





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

    % loop through potential list of files
    armaData = struct;
    lengthData = struct;

    for iFile = 1:length(strainFile)

        dataSuffix = regexp(strainFile{iFile},'strains_([\w]+).mat','tokens');
        dataSuffix = char(dataSuffix{1});
        % armaData is a collection of many files, while strainInfo is already
        % the correct structure that contains the optimal model for each
        armaDataTmp = load(fullfile(dataPath,sprintf('resFitVelAndLen_%s',dataSuffix)));

        % this would be really annoying without a loop
        fn = fieldnames(armaDataTmp);
        for i=1:length(fn)
            armaData.(fn{i}) = armaDataTmp.(fn{i});
        end

        % try load length file
        try
            lengthDataTmp = load(fullfile(dataPath,sprintf('lengthSeries_%s',dataSuffix)));
        catch
            lengthDataTmp = struct;
        end
        fn = fieldnames(lengthDataTmp);
        for i=1:length(fn)
            lengthData.(fn{i}) = lengthDataTmp.(fn{i});
        end


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
        goodIdx = find(goodIdx);
        for i=1:length(goodIdx)
            data(i).lengthSeries = ...
                lengthData.(sprintf('length%s',strainInfo(selectionIdx(goodIdx(i))).name));
            if options.multiply
                data(i).lengthSeries = repmat(data(i).lengthSeries,[1,options.multiply]);
                data(i).numObserve = data(i).numObserve * options.multiply;
            end
        end
    end

end % load data



%============================


%============================
%% Group according to WNV
%============================

groupingOptions.labels = strvcat(data.name);
distanceFunctionParameters.cutoff = options.wnv1_cutoff;
distanceFunctionParameters.mode = options.wnv1_mode;
[links, groupIdx, groupedData, linkData] = ...
    groupData(data,'groupArma_distance_WNV',...
    groupingOptions,distanceFunctionParameters);


% label axis
set(get(groupedData.plotHandles.axesH,'XLabel'),'String','-^1^0log(probability)')
set(groupedData.plotHandles.figureH,'Name','WNV - Round 1');

% have the tree grow from the right
set(groupedData.plotHandles.axesH,'XDir','reverse','YAxisLocation','right')


%===========================


%===========================
%% Group according to ARMA
%===========================

%options.mode = 2;
distanceFunctionParameters.cutoff = options.arma_cutoff;
distanceFunctionParameters.mode = options.arma_mode;

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
        groupedData.collectedData(iGroup).data,'groupArma_distance_ARMA',...
        groupingOptions,distanceFunctionParameters);

    % steal all we can from the figure
    ah = groupedData.collectedData(iGroup).subGroups.plotHandles.axesH;
    % steal lines
    copyobj(get(ah,'Children'),sh)
    % steal y-labels
    set(sh,'YTick',get(ah,'YTick'),'YTickLabel',get(ah,'YTickLabel'),'YLim',get(ah,'YLim'))
    set(get(sh,'XLabel'),'String','-^1^0log(probability)')
    % grow tree from the right
    set(sh,'XDir','reverse','YAxisLocation','right')

    % close figure
    delete(groupedData.collectedData(iGroup).subGroups.plotHandles.figureH);

end

%===========================
%% Group according to WNV
%===========================

%options.mode = 2;
distanceFunctionParameters.cutoff = options.wnv2_cutoff;
distanceFunctionParameters.mode = options.wnv2_mode;


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
                groupedData.collectedData(iGroup).subGroups.collectedData(jGroup).data,...
                'groupArma_distance_WNV',groupingOptions,distanceFunctionParameters);
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

            % grow tree from the right
            set(sh,'XDir','reverse','YAxisLocation','right')
            % close figure
            delete(groupedData.collectedData(iGroup).subGroups.collectedData(jGroup).subGroups.plotHandles.figureH);
        end
    end

end


%=============================
%% CREATE OUTPUT
%=============================

% the whole groupedData-thing is a tad complicated. Return a
% groupIndex-list instead that has as many rows as data sets, and that has
% a column for every grouping step that indicates to which group the data
% belongs

% counters of group-numbers already used
jLevelCt = 0;
kLevelCt = 0;

groupIndex = zeros(nData,3);

for iGroup = 1:length(groupedData.collectedData)

    % check whether there are any subgroups
    if isempty(groupedData.collectedData(iGroup).subGroups)
        % all falls into the same group
        jGroup = 1;
        kGroup = 1;
        % update kLevelCt here, too!
        kLevelCt = kLevelCt + kGroup;

        % fill i,j,k
        groupIndex(groupedData.collectedData(iGroup).dataIdx,:) = ...
            repmat([iGroup, jGroup+jLevelCt, kGroup+kLevelCt],...
            length(groupedData.collectedData(iGroup).dataIdx),1);

    else % there are subgroups

        for jGroup = 1:length(groupedData.collectedData(iGroup).subGroups.collectedData)

            % check whether there are any subgroups
            if isempty(groupedData.collectedData(iGroup).subGroups.collectedData(jGroup).subGroups)
                % all falls into the same group
                kGroup = 1;

                % fill i,j,k
                dataIdx = groupedData.collectedData(iGroup).dataIdx(...
                    groupedData.collectedData(iGroup).subGroups.collectedData(jGroup).dataIdx);
                groupIndex(dataIdx,:) = ...
                    repmat([iGroup, jGroup+jLevelCt, kGroup+kLevelCt],...
                    length(dataIdx),1);

            else % there are subgroups

                for kGroup = 1:length(groupedData.collectedData(iGroup).subGroups.collectedData(jGroup).subGroups.collectedData)
                    
                    % fill i,j,k
                    dataIdx = groupedData.collectedData(iGroup).dataIdx(...
                    groupedData.collectedData(iGroup).subGroups.collectedData(jGroup).dataIdx(...
                    groupedData.collectedData(iGroup).subGroups.collectedData(jGroup).subGroups.collectedData(kGroup).dataIdx));
                groupIndex(dataIdx,:) = ...
                    repmat([iGroup, jGroup+jLevelCt, kGroup+kLevelCt],...
                    length(dataIdx),1);
                    
                end

            end

            % update number of 2nd-level groups
            kLevelCt = kLevelCt + kGroup;

        end % loop 2nd level


    end

    % update number of 2nd-level groups
    jLevelCt = jLevelCt + jGroup;

end % loop 1st level=======

%=============================
%% CLEANUP
setappdata(0,'UseNativeSystemDialogs',sysDialogState);
