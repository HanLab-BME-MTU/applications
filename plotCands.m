function plotCands(candsCoords,candsPerFrame,first,last)

% get dimensions of roi where cands data comes from
[name path]=uigetfile('.tif','Select one of the roi images used to get cands');
imInfo = imfinfo([path name]);
roiWidth = imInfo.Width
roiHeight = imInfo.Height

numFiles = length(candsPerFrame);
if first<=0 | last>numFiles
    disp('first and/or last frame is out of range')
end
if first>last
    disp('last frame must be greater than first')
end

% make origin in upper left corner to match images
set(gca,'YDir','reverse');
for i=first:last
    scatter(candsCoords(1:candsPerFrame(i),2,i),candsCoords(1:candsPerFrame(i),1,i),3,[rand rand rand]);
    hold on
end
axes('XAxisLocation',top,'XLim',[0 roiWidth],'YLim',[0,roiHeight])
hold off
