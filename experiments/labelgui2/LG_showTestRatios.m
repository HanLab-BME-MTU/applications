function figureHandle = LG_showTestRatios(testRatios, dataProperties,figurePosition)
%LG_showTestRatios is the routine to display testRatios and the current cutoff
% figurePosition is optional

% read amplitude cutoff
cutValue = dataProperties.amplitudeCutoff;

% create figure
figureHandle = figure('Name',sprintf('cutoff for %s',dataProperties.name));
% place correctly if there is a remembered position
if nargin > 2 && ~isempty(figurePosition)
    set(figureHandle,'Position',figurePosition);
end

% read testRatios
tmp = cat(1,testRatios{:});
times = (tmp(:,1));
ratios = tmp(:,2);
nTimepoints = max(times);

% plot testRatios
t1idx = ismember(times,1:3:nTimepoints);
t2idx = ismember(times,2:3:nTimepoints);
t3idx = ismember(times,3:3:nTimepoints);
h1=plot(times(t1idx),ratios(t1idx),'+r',...
    times(t2idx),ratios(t2idx),'+g',...
    times(t3idx),ratios(t3idx),'+b');

% plot cutoff
hold on
h2=plot([1,nTimepoints],[cutValue,cutValue],'r');
  
% allow navigating to individual frames by click
set([h1;h2;gca],'ButtonDownFcn','LG_gotoFrameWithClick');