function [groupIndex,groupedData] = groupArmaDescriptors(data,options,verbose)
%GROUPARMADESCRIPTORS groups models according to similarity of their ARMA coefficients
%
% SYNOPSIS: [groupIndex,groupedData] = groupArmaDescriptors(data,options,verbose)
%
% INPUT     data (opt): structure similar to the output of armaxFitKalman, with
%                 three additional fields: 
%                   - orderLen OR orderVel (number of fit parameters), 
%                   - type: 'Len' or 'Vel' depending on whether length or
%                     velocities were fitted.
%                   - name (name of the dataSet). 
%                 If empty, the code will ask for a file called
%                 strain_???.mat, an then look for a file called
%                 resFitVelAndLen_???.mat, and, possibly,
%                 lengthSeries_???.mat OR velocitySeries_???.mat
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
%                 type         'Len' or 'Vel' - whether to group length or
%                              velocity series. Default: 'Len'
%           verbose : [0/{1}] whether or not to display figures
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
    groupIndex = [];
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
def_options.wnv1_cutoff = -2;
def_options.wnv1_mode = [0 1e-12];
def_options.arma_cutoff = 5e-4;
def_options.arma_mode = [0,5e-4];
def_options.wnv2_cutoff = 1e-12;
def_options.wnv2_mode = [0 1e-12];
def_options.plot = 1; % plot results
def_options.type = 'Len'; % length series as default
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

% check case of type
options.type = lower(options.type);
options.type(1) = upper(options.type(1));


% check verbose
if nargin < 3 || isempty(verbose)
    verbose = 1;
end

if nargin == 0 || isempty(data)
    data = groupArma_loadData(options);
end % load data

%get number of data sets
nData = length(data);

if nData == 0
    disp('--no data loaded. groupArmaDescriptors aborted')
    return
end
    


%============================


%============================
%% Group according to WNV
%============================

groupingOptions.labels = strvcat(data.name);
groupingOptions.verbose = verbose;
distanceFunctionParameters.cutoff = options.wnv1_cutoff;
distanceFunctionParameters.mode = options.wnv1_mode;
[links, groupIdx, groupedData, linkData] = ...
    groupData(data,'groupArma_distance_WNV',...
    groupingOptions,distanceFunctionParameters);

if verbose
    % label axis
    set(get(groupedData.plotHandles.axesH,'XLabel'),'String','-^1^0log(probability)')
    set(groupedData.plotHandles.figureH,'Name','WNV - Round 1');

    % have the tree grow from the right
    set(groupedData.plotHandles.axesH,'XDir','reverse','YAxisLocation','right')
end

%===========================


%===========================
%% Group according to ARMA
%===========================

%options.mode = 2;
distanceFunctionParameters.cutoff = options.arma_cutoff;
distanceFunctionParameters.mode = options.arma_mode;

nGroups = length(groupedData.collectedData);

% plot
if verbose
    plotRows = ceil(sqrt(nGroups));
    plotCols = ceil(nGroups/plotRows);
    fh = figure('Name','ARMA-clustering');
end

for iGroup = 1:nGroups

    % group, find subgroups
    groupingOptions.labels = strvcat(groupedData.collectedData(iGroup).data.name);
    groupingOptions.verbose = verbose;

    [groupedData.collectedData(iGroup).links,groupedData.collectedData(iGroup).groupIdx,...
        groupedData.collectedData(iGroup).subGroups,groupedData.collectedData(iGroup).linkData] = groupData(...
        groupedData.collectedData(iGroup).data,'groupArma_distance_ARMA',...
        groupingOptions,distanceFunctionParameters);

    if verbose
        figure(fh);
        sh = subplot(plotRows,plotCols,iGroup);
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
end

%===========================
%% Group according to WNV
%===========================

%options.mode = 2;
distanceFunctionParameters.cutoff = options.wnv2_cutoff;
distanceFunctionParameters.mode = options.wnv2_mode;


for iGroup = 1:nGroups

    nSubGroups = length(groupedData.collectedData(iGroup).subGroups.collectedData);

    if verbose
        % plot
        plotRows = ceil(sqrt(nSubGroups));
        plotCols = ceil(nSubGroups/plotRows);
        fh = figure('Name',sprintf('WNV - Round 2. Arma-Group %i',iGroup));
    end

    for jGroup = 1:nSubGroups

        if verbose
            figure(fh);
            sh = subplot(plotRows,plotCols,jGroup);
        end

        % group
        if length(groupedData.collectedData(iGroup).subGroups.collectedData(jGroup).data)>1
            groupingOptions.labels = strvcat(groupedData.collectedData(iGroup).subGroups.collectedData(jGroup).data.name);
            groupingOptions.verbose = verbose;

            [groupedData.collectedData(iGroup).subGroups.collectedData(jGroup).links,...
                groupedData.collectedData(iGroup).subGroups.collectedData(jGroup).groupIdx,...
                groupedData.collectedData(iGroup).subGroups.collectedData(jGroup).subGroups,...
                groupedData.collectedData(iGroup).subGroups.collectedData(jGroup).linkData] = groupData(...
                groupedData.collectedData(iGroup).subGroups.collectedData(jGroup).data,...
                'groupArma_distance_WNV',groupingOptions,distanceFunctionParameters);
            if verbose
                xlabel('-^1^0log(probability)')
            end
        end

        if isfield(groupedData.collectedData(iGroup).subGroups.collectedData(jGroup),'subGroups') &&...
                ~isempty(groupedData.collectedData(iGroup).subGroups.collectedData(jGroup).subGroups) && verbose
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
    if ~isfield(groupedData.collectedData(iGroup),'subGroups') || ...
            isempty(groupedData.collectedData(iGroup).subGroups)
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
            if ~isfield(groupedData.collectedData(iGroup).subGroups.collectedData(jGroup),'subGroups') || ...
                    isempty(groupedData.collectedData(iGroup).subGroups.collectedData(jGroup).subGroups)
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
