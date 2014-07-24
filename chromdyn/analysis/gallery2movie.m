function gallery2movie(fh)
%GALLERY2MOVIE transforms the spindle projection galleries to .avi movies
%
% SYNOPSIS: gallery2movie(fh)
%
% INPUT fh: handle to figure with the gallery
%
% OUTPUT 
%
% REMARKS
%
% SEE ALSO
%

% created with MATLAB ver.: 7.11.0.584 (R2010b) on Mac OS X  Version: 10.6.4 Build: 10F569 
%
% created by: Jonas Dorn
% DATE: 04-Oct-2010
%
% Last revision $Rev: 1507 $ $Date: 2010-10-06 09:53:26 -0400 (Wed, 06 Oct 2010) $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% we need to
% 1. get cData, xData and yData from the images
% 2. get axes limits from the corresponding axes
% 3. init movie file with videoWriter
% 4. Loop to create figures and call getframe on the axes

if nargin < 1 || isempty(fh)
    fh = gcf;
end

% get images
imageH = findall(fh,'type','image');

nImages = length(imageH);

plotData(1:nImages) = struct('cData',[],'xData',[],'yData',[],'spindleLength',[]);

for im = 1:nImages
    plotData(im).cData = get(imageH(im),'cData');
    plotData(im).xData = get(imageH(im),'xData');
    plotData(im).yData = get(imageH(im),'yData');
    titleStr = get(get(get(imageH(im),'parent'),'title'),'string');
    plotData(im).spindleLength = str2double(strtok(titleStr{1},'--'));
end

% get axes limits
xlim = get(get(imageH(1),'Parent'),'xlim');
ylim = get(get(imageH(1),'Parent'),'ylim');

% get title
fileName = titleStr{2}(1:end-4);

% create videoWriter object
vObj = VideoWriter(fileName);
vObj.FrameRate = 3;
vObj.Quality = 75;

% open for writing
vObj.open;

% loop and plot. Count backwards because handles to longest spindles are
% first (since they were plotted last)
for im = nImages:-1:1
    fhx = figure;
    
    ih = imshow(plotData(im).cData,...
        'XData',plotData(im).xData,'YData',plotData(im).yData);
    set(get(ih,'parent'),'xlim',xlim,'ylim',ylim,'visible','on','color','k','position',[0.1,0.1,0.8,0.8])
    fullscreen(fhx)
    % add spindle length
    text(xlim(1)*0.95,ylim(2)*0.90,num2str(plotData(im).spindleLength),...
        'color','w','fontSize',48);
    
    drawnow
    
    currFrame = getframe(get(ih,'parent'));
    vObj.writeVideo(currFrame);
    
    close(fhx)
end

vObj.close;

