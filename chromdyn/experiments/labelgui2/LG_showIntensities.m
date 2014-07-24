function [figureHandle,objectHandles] = LG_showIntensities(idlist,idlistData,colorMap,figurePosition)
%LG_SHOWINTENSITIES shows a plot of spot and tag intensities in the chromdyn project
%
% SYNOPSIS: figureHandle = LG_showIntensities(idlist,idlistData,figurePosition)
%
% INPUT     idlist - idlist
%           idlistData - (opt) idlistData generated with LG_readIdlist
%                              if no idlistData, supply dataProperties here
%                              instead
%           colorMap - (opt) colorOrder of the tags
%           figurePosition - (opt) requested figure position (if it's a
%                   figure handle, plot will be into existing axes)
%
% OUTPUT    figureHandle - handle to the intensity figure
%
% REMARKS The function will create one set of axes to show the intensity
%         fit, and one set of axes to show the tag- and spot intensities
%
% created with MATLAB ver.: 7.1.0.246 (R14) Service Pack 3 on Windows_NT
%
% created by: jdorn
% DATE: 28-Mar-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%================
% TEST INPUT
%================

if nargin < 1 || isempty(idlist)
    error('not enough input arguments')
end

% read idlistdata if necessary. It will calculate lots of unnecessary data,
% but it's reasonably fast
if nargin < 2 || isempty(idlistData) 
    error('please specify idlistData or dataProperties as second input argument')
elseif ~isfield(idlistData,'nSpots')
    idlistData = LG_readIdlistData(idlist, idlistData);
end

% colormap: allow empty
if nargin < 3
    colorMap = [];
end

% figurePosition: allow empty
if nargin < 4
    figurePosition = [];
end

%========================


%==================
% LAUNCH FIGURE
%==================

% only launch figure if necessary
if isempty(figurePosition) || length(figurePosition) > 1 || ~ishandle(figurePosition)
    figureHandle = figure('Name',sprintf('Intensities for %s',idlist(1).stats.name));
    % place correctly if there is a remembered position
    if  ~isempty(figurePosition)
        set(figureHandle,'Position',figurePosition);
    end
else
    figureHandle = figurePosition;
end
figure(figureHandle);

%==================


%=================
% SHOW INT-FIT
%=================

% show intFit on top
ah(1,1) = subplot(2,1,1,'replace');

% get data
time = 1:length(idlist);
fitParms = idlist(1).stats.intFit.xFit;
% we used to have first fitParms as log, therefore ensure backwards
% compatibility
if any(fitParms(1:end-1) < 0)
yFit = exp(fitParms(1:end-1)) * exp(fitParms(end) * time);
else
    yFit = (fitParms(1:end-1)) * exp(fitParms(end) * time);
end
yFit = yFit';

nSpotList = unique(idlistData.nSpots);
nSpotList(nSpotList == 0) = [];


% store handles
hfit = zeros(length(nSpotList)+length(fitParms)-1,1);

% plot the fitted lines
hfit(1:length(fitParms)-1)=plot(ah(1),time,yFit);
colorOrder = get(ah(1),'ColorOrder');
hold on

% plot the dots in the correct color
for i = 1:length(nSpotList)
    plotIdx = find(idlistData.nSpots == nSpotList(i));
    hfit(i+length(fitParms)-1) = plot(ah(1),time(plotIdx),idlist(1).stats.ampList(plotIdx),...
        '.','Color',colorOrder(wraparound(i,[1;size(colorOrder,1)]),:));

end

% label the plot, add legend
legendString = [num2str(nSpotList),repmat(' spots',length(nSpotList),1)];
legend(ah(1),legendString);
set(get(ah(1),'Title'),'Interpreter','none','String',...
    sprintf('time constant of intensity fit %1.5f',fitParms(end)));

%=================

%======================
% SHOW TAG INTENSITIES
%======================

% check for good/"bad" tags
goodTagIdx = find(idlist(idlistData.goodTimes(1)).linklist(:,5) < 2);
badTagIdx = find(idlist(idlistData.goodTimes(1)).linklist(:,5) > 1);

% make second axes, make colors in case we don't get them passed down.
ah(2,1) = subplot(2,1,2,'replace');
if isempty(colorMap)
    colorMap = hsv(idlistData.maxTags);
end
set(ah(2),'NextPlot','add');
hold on
% plot good tags
% store handles (there'll be 3 handles per tag!)
hg = zeros(3*length(goodTagIdx),1);
for iTag = 1:length(goodTagIdx)
    % read list of intensities, amplitudes, isSpot
    currentTagIdx = goodTagIdx(iTag);
    intList = catStruct(1,sprintf('idlist.linklist(%i,[1,8,3])',currentTagIdx));
    isSpot = intList(:,3) ~= 1;

    % plot the entire course of intensities first, then distinguish between
    % true and estimated spot via markers
    hg((iTag-1)*3+1) = plot(ah(2),intList(:,1),intList(:,2),'-');
    hg((iTag-1)*3+2) = plot(ah(2),intList(isSpot,1),intList(isSpot,2),'o');
    % if there's nothing to plot, plot returns []!
    if any(~isSpot)
        hg((iTag-1)*3+3) = plot(ah(2),intList(~isSpot,1),intList(~isSpot,2),'*');
    else
        hg((iTag-1)*3+3) = plot(ah(2),0,0,'-');
    end
    
    % make sure the color is correct
    set(hg((iTag-1)*3+1:iTag*3),'Color',colorMap(currentTagIdx,:));
end

% plot estimated tags
hb = zeros(length(badTagIdx),1);
for iTag = 1:length(badTagIdx)
    % read list of intensities, amplitudes, isSpot
    currentTagIdx = badTagIdx(iTag);
    intList = catStruct(1,sprintf('idlist.linklist(%i,[1,8,3])',currentTagIdx));
    isSpot = intList(:,3) ~= 1;

    % plot only the individual occurence
    hb(iTag) = plot(ah(2),intList(isSpot,1),intList(isSpot,2),'h');
    % make sure the color is correct
    set(hb(iTag),'Color',colorMap(currentTagIdx,:));
end

% name axes
set(get(ah(2),'Title'),'Interpreter','none','String','tag intensities');

% store handles to make everything clickable
objectHandles = [hfit;hg;hb;ah];
%======================

