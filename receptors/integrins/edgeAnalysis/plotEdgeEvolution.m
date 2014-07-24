function plotEdgeEvolution(MD)
%PLOTEDGEEVOLUTION makes a plot of all edges in a movie with time color-coding
%
%SYNPOSIS plotEdgeEvolution(MD)
%
%INPUT  MD          : movieData with cell mask.
%
%OUTPUT the plot.
%
%Khuloud Jaqaman, October 2013

%% preamble

%get images
imageDir = MD.channels_.channelPath_;
imageFileListing = dir([imageDir filesep '*.tif']);
if isempty(imageFileListing)
    imageFileListing = dir([imageDir filesep '*.tiff']);
end

%get masks
maskDir = MD.processes_{2}.outFilePaths_{1};
maskFileListing = dir([maskDir filesep '*.tif']);
if isempty(maskFileListing)
    maskFileListing = dir([maskDir filesep '*.tiff']);
end

%get number of frames
numFramesMovie = length(imageFileListing);

%% make plot

%generate color-coding for time
if mod(numFramesMovie,2) == 1
    midPoint = ceil(numFramesMovie/2);
    colorUp = linspace(0,1,midPoint)';
    colorDn = colorUp(end:-1:1);
    edgeColor = [zeros(midPoint,1) colorDn colorUp; colorUp(2:end) zeros(midPoint-1,1) colorDn(2:end)];
else
    midPoint = numFramesMovie/2;
    colorUp = linspace(0,1,midPoint)';
    colorDn = colorUp(end:-1:1);
    edgeColor = [zeros(midPoint,1) colorDn colorUp; colorUp zeros(midPoint,1)];
end

%display first image
imageStack = imread(fullfile(imageDir,imageFileListing(1).name));
figure
imshow(imageStack,[prctile(imageStack(:),5) prctile(imageStack(:),95)]);
% hold on

%go over all frames and plot mask boundaries
figure
imshow(ones(size(imageStack)),[]);
hold on
for iFrame = 1 : numFramesMovie
    maskStack1 = imread(fullfile(maskDir,maskFileListing(iFrame).name));
    maskBounds = bwboundaries(maskStack1);
    cellfun(@(x)(plot(x(:,2),x(:,1),'Color',edgeColor(iFrame,:),'LineWidth',2)),maskBounds);    
end


%% ~~~ end ~~~

