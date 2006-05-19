function [figureHandle, objectHandles] = LG_showDistances(idlist, dataProperties, figurePosition)
%LG_SHOWDISTANCES shows distance plots between up to four tags
%
% SYNOPSIS: [figureHandle, objectHandles] = LG_showDistances(idlist, dataProperties)
%
% INPUT   idlist - new idlist
%         dataProperties - dataProperties structure (see
%                          defaultDataProperties for more information)  
%         figurePosition (opt) - figure position
%
% OUTPUT figureHandle: handle to the figure
%        objectHandles: handles to axis, plots
%
% REMARKS
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 27-Apr-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=================================
%% TEST INPUT AND GET BASIC STUFF
%=================================

if nargin < 2 || isempty(idlist) || isempty(dataProperties)
    error('LG_showDistances needs at least two input arguments');
end

if nargin < 3 
    figurePosition = [];
end

% get idlist info
nTimepoints = length(idlist);

% calculate distances
[distances,dummy,dummy,dummy,idlist,idxLists] = idlist2distMat(idlist, dataProperties,0,1);

% check for number of tags. We can display up to 5, above that, user has to
% select
tagList = idlist(1).stats.labelcolor;
nTags = length(tagList);

if nTags > 5
    h = helpdlg('please select up to 5 tags between which the distance will be plotted');
    uiwait(h);
    [tagIdx, tagList] = listSelectGUI(idlist(1).stats.labelcolor,5);
    nTags = length(tagIdx);
else
    tagIdx = 1:nTags;
end

% select plot distribution
switch nTags
    case {0,1}
        error('no tags for plotting!')
    case {2,3,4,5}
        subAxes = [max(nTags-2,1),nTags-1];
    otherwise
        error('too many tags!')
end

% make pair index. 
[u,v] = ndgrid(1:nTags,1:nTags);
idx = find(tril(u,-1));
pairIdx = [u(idx) v(idx)];
nPairs = length(idx);

%===============================


%===============================
%% PLOT
%===============================

figureHandle = figure('Name',sprintf('Distances for %s',idlist(1).stats.name));
% place correctly if there is a remembered position
if  ~isempty(figurePosition)
    set(figureHandle,'Position',figurePosition);
end

% collect distances here so that we can get the maximum
dist2plot = zeros(nTimepoints,nPairs);
for i=1:nPairs
    dist2plot(:,i) = squeeze(distances(tagIdx(pairIdx(i,1)),...
        tagIdx(pairIdx(i,2)),:));
    estimatedIdx(:,i) = any(idxLists.estimatedTag(:,tagIdx(pairIdx(i,:))),2);
end
maxDist = nanmax(dist2plot(~estimatedIdx(:)));

% loop through tags and plot
objectHandles = zeros(4*nPairs);

for i=1:nPairs
    % make axes
    ah = subplot(subAxes(1),subAxes(2),i);
    objectHandles((i-1)*4+1) = ah;
    % plot distance
    x = (1:nTimepoints)';
    y = dist2plot(:,i);
    % plot a broken line everywhere, then add a solid line for good
    % distances.
    xGood = x;
    xGood(estimatedIdx(:,i)) = NaN;
    yGood = y;
    yGood(estimatedIdx(:,i)) = NaN;
    objectHandles((i-1)*4+2) = plot(x,y,'--');
    hold on
    objectHandles((i-1)*4+3) = plot(xGood, yGood,'-o');
    objectHandles((i-1)*4+4) = plot(x(estimatedIdx(:,i)),y(estimatedIdx(:,i)),'*');
    
    % write title
    set(get(ah,'Title'),'Interpreter','none','String',...
        sprintf('%s - %s',tagList{tagIdx(pairIdx(i,1))},tagList{tagIdx(pairIdx(i,2))}));
    ylim([0,maxDist]);
end