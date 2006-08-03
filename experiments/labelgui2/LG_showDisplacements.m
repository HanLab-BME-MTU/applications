function [figureHandle, objectHandles] = LG_showDisplacements(idlist, dataProperties, figurePosition, colorMap)
%LG_SHOWDISPLACEMENTS shows displacement of up to 25 tags
%
% SYNOPSIS: [figureHandle, objectHandles] = LG_showDisplacements(idlist, dataProperties)
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

% colormap: because of possible tag selection, we can't assign default here
if nargin < 4 
    colorMap = [];
end

% get idlist info
nTimepoints = length(idlist);

% calculate distances
[dummy,dummy,displacements,dummy,idlist] = idlist2distMat(idlist, dataProperties);
% don't remember sigma
displacements = displacements(:,:,1);

% check for number of tags. We can display up to 5, above that, user has to
% select
tagList = idlist(1).stats.labelcolor;
nTags = length(tagList);

if nTags > 25
    h = helpdlg('please select up to 25 tags between which the distance will be plotted');
    uiwait(h);
    [tagIdx, tagList] = listSelectGUI(idlist(1).stats.labelcolor,25);
    nTags = length(tagIdx);
    if ~isempty(colorMap)
        colorMap = colorMap(tagIdx,:);
    end
else
    tagIdx = 1:nTags;
end

% select plot distribution
if nTags == 0
    error('not enough tags to plot!')
elseif nTags > 25
    error('too many tags to plot!')
else
    nRows = floor(sqrt(nTags));
    nCols = ceil(nTags/nRows);
end


%===============================


%===============================
%% PLOT
%===============================

figureHandle = figure('Name',sprintf('Displacements for %s',idlist(1).stats.name));
% place correctly if there is a remembered position
if  ~isempty(figurePosition)
    set(figureHandle,'Position',figurePosition);
end

% collect distances here so that we can get the maximum
disp2plot = displacements(:,tagIdx);

maxDisp = nanmax(disp2plot(:));

% loop through tags and plot
objectHandles = zeros(2*nTags);

% get colormap
if isempty(colorMap)
    colorMap = repmat([0,0,1], nTags, 1);
end

for i=1:nTags
    % make axes
    ah = subplot(nRows,nCols,i);
    objectHandles((i-1)*2+1) = ah;
    % plot distance
    objectHandles(i*2) = plot(disp2plot(:,i),'d-','Color',colorMap(i,:));
    % write title
    set(get(ah,'Title'),'Interpreter','none','String',...
        sprintf('%s',tagList{(i)}));
    ylim([0,maxDisp]);
end