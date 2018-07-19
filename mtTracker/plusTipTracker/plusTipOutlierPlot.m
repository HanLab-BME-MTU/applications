function [groupListFilt,outlierIdx,cMap]=plusTipOutlierPlot(plusTipDataset,groupList,paramName,doPlot,cMap,doRem)
% plusTipOutlierPlot removes outliers from a list of projects and plots their values
%
% SYNOPSIS: 
% [groupListFilt,outlierIdx,cMap]=plusTipOutlierPlot(plusTipDataset,groupList,paramName,doPlot,cMap,doRem)
%
% INPUT
%   plusTipDataset : dataset array from plusTipGetGroupDataset, which
%                    contains at least (but may contain more than) all the
%                    movies in groupList
%   groupList      : n x 2 cell array containing group labels in the first
%                    column and project paths in the second. this will
%                    usually be the output of plusTipPickGroups
%   paramName      : parameter name from plusTipDataset for which outlier
%                    detection and plotting should be done. if empty, user
%                    will be asked to select from a list.
%   doPlot         : (1) to plot the parameter values by group, 0 to skip
%   cMap (opt)     : nGroups x 3 color map, if not given random colors will
%                    be chosen. note that the chosen cMap is also an output
%                    in case you want to plot data from a different
%                    parameters using the same colors for the groups
%   doRem          : (1) to determine and remove outliers, 0 to skip.
%
% OUTPUT
%   groupListFilt  : same as groupList, except the outlier projects are
%                    removed
%   outlierIdx     : index list corresponding to the input groupList where
%                    outliers were found
%   cMap           : color map used, see input explanation above
%
%   if doPlot=1, a figure will be greated showing the parameter value for
%   each group. the groups will be distinguished by color and shape.
%   outliers are marked by a red square.


homeDir=pwd;
if nargin<1 || isempty(plusTipDataset)
    [fileName,pathName] = uigetfile('*.mat','Select plusTipDataset.mat file');
    if fileName==0
        display('no file selected')
        return
    end
    load([pathName filesep fileName]);
    cd(pathName)
end

if nargin<2 || isempty(groupList)
    groupList=combineGroupListFiles;
end
cd(homeDir)

if nargin<3 || isempty(paramName)
    % let the user pick the parameter
    paramNameList=plusTipDataset.Properties.VarNames;
    paramNameList(1:13)=[]; % these come from tracking parameters
    selection=listSelectGUI(paramNameList,1,'copy',1);
    paramName=paramNameList{selection};
end

if nargin<4 || isempty(doPlot)
    doPlot=1;
end

if nargin<5 || isempty(cMap)
    cMap=[];
end

if nargin<6 || isempty(doRem)
    doRem=1;
end


groupListFilt=groupList;

projGroupName=groupList(:,1);
projGroupDir=cellfun(@(x) formatPath(x),groupList(:,2),'uniformoutput',0);


% fix the names if there are spaces or hyphens and append prefix 'grp'
projGroupName=cellfun(@(x) strrep(x,'-','_'),projGroupName,'uniformoutput',0);
projGroupName=cellfun(@(x) strrep(x,' ','_'),projGroupName,'uniformoutput',0);
projGroupName=cellfun(@(x) ['grp_' x],projGroupName,'uniformoutput',0);


% count unique groups and keep them in order of the original list
[btwGrpNames,m,projGroupIdx] = unique(projGroupName);
[b,idx]=sort(m);
btwGrpNames=btwGrpNames(idx);
nGroups=length(btwGrpNames);


outlierIdx=[];
c=1;
gsAll=zeros(size(groupList,1),1);

for iGroup = 1:nGroups

    % indices of projects in iGroup
    tempIdx=strmatch(btwGrpNames(iGroup),projGroupName,'exact');

    % find the data index in plusTipDataset corresponding to each project
    % in the iGroup
    datasetIdx=zeros(length(tempIdx),1);
    for iProj=1:length(tempIdx)
        datasetIdx(iProj,1)=find(cellfun(@(x) ~isempty(strmatch(formatPath(x),formatPath(projGroupDir{tempIdx(iProj)}),'exact')),plusTipDataset.projPath));
    end

    % get the parameter values
    paramValues=plusTipDataset.(paramName)(datasetIdx);

    if doPlot==1
        if iGroup==1
            % properties for plotting
            nRep=ceil(nGroups/3);
            mrkTpe=repmat(['o';'^';'s'],nRep,1); % we'll make this into a cell later
            prop_name(1) = {'Marker'};
            prop_name(2) = {'MarkerFaceColor'};
            prop_name(3) = {'MarkerEdgeColor'};

            % set the colormap to be random
            if isempty(cMap)
                cM=hsv(nGroups);
                cM=cM(randsample(nGroups,nGroups),:);
                cMap=mat2cell(cM,ones(nGroups,1),3);
            end

            figure; hold on
        end

        nC=length(paramValues);
        prop_values(1:nC,1) = {mrkTpe(iGroup)};
        prop_values(1:nC,2) = cMap(iGroup,:);
        prop_values(1:nC,3) = cMap(iGroup,:);

        % use plot instead of scatter so more flexibility with properties. to
        % do this, make 2 x nPoints matrix where the second row is all NaNs and
        % then use the plot function
        xCoord=[tempIdx nan(size(tempIdx))]';
        yCoord=[paramValues nan(size(paramValues))]';
        h=plot(xCoord,yCoord,'.');
        set(h,prop_name,prop_values)

        % store the handle for the first member of each group for the legend
        h1(iGroup)=h(1);

        % these have to be defined each time
        clear prop_values

        % keep all the parameter values in a vector
        gsAll(c:c+length(paramValues)-1,1)=paramValues;

        c=c+length(paramValues);
    end
    if doRem==1
        % find outliers and plot them
        outliers=find(rousseeuwCrit(paramValues));
        if ~isempty(outliers)
            outlierIdx=[outlierIdx; tempIdx(outliers)];
        end
    end



end
if doPlot==1
    scatter(outlierIdx,gsAll(outlierIdx),[],'rs','SizeData',(72/8)^2)
    btwGrpNames=cellfun(@(x) strrep(x,'_','-'),btwGrpNames,'uniformoutput',0);
    legend(h1,btwGrpNames,'location','BestOutside');
    xlabel('groupList index')

    str=strrep(paramName,'_std','');
    str=strrep(str,'_',' ');

    ylabel(str)
end
if doRem==1
    groupListFilt(outlierIdx,:)=[];
end