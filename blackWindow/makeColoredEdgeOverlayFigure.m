function figHan = makeColoredEdgeOverlayFigure(MD,iFrame,iChan,colorMap)
%MAKECOLOREDEDGEOVERLAYFIGURE makes the figure with color-coded edge timelapse overlaid on a single image
% 
% makeColoredEdgeOverlayFigure;
% makeColoredEdgeOverlayFigure(movieData,frameNumber,channelNumber,@colorMap);
% 
% Displays one image from the input movie with the edge outlines from every
% frame of the movie overlaid on this image, color-coded for time.
% 
% Input: (All inputs are optional. When necessary the user will be asked.)
%
%       movieData : MovieData object for the movie to make figure from.
%                   Optional. If not input, user will be asked to select
%                   a movieData file.
%
%      frameNubmer: Frame number to display image from.
%                   Optional. If not input, user will be asked.
%
%    channelNumber: Channel number to display image from.
%                   Optional. If not input, user will be asked.
%
%       @colorMap : Handle to colormap to use for display, e.g. @jet or @hot
%                   Optional. Default is to use @jet.
%
%
% Hunter Elliott
% 9/2012
% 


%% ------------ Input ----------- %%

if nargin < 1 || isempty(MD);
    [fileName,filePath] = uigetfile('*.mat','Please select a movieData file:');        
    if fileName == 0
        return
    else
        MD = MovieData.load([filePath filesep fileName]);    
    end
else
    MD.sanityCheck;
end
   
iProt = MD.getProcessIndex('ProtrusionProcess',1,1);
if isempty(iProt)
    error('This movie does not have a valid protrusion process! Please run protrusion calculation first!')
end

if nargin < 2 || isempty(iFrame);    
    iFrame = str2double(inputdlg('Input a frame number to overlay edges on:','Frame Number'));
end

nFrames = MD.nFrames_;
if isempty(iFrame) || iFrame > nFrames || iFrame <= 0
    error(['You must select a valid frame number to confinue!! Frame number must be between 1 and ' num2str(MD.nFrames_) '!!'])
end

nChan = numel(MD.channels_);
if nChan == 1
    iChan = 1;
    OK = true;
elseif nargin < 3 || isempty(iChan)
    [iChan OK] = listdlg('ListString',cellfun(@num2str,1:nChan),'SelectionMode','single',...
        'PromptString','Select a channel to overlay edges on:');
end
if ~OK
    return
end

if nargin < 4 || isempty(colorMap)
    colorMap = @jet;
end


%% ----------------- Init ---------- %%

imNames = MD.getImageFileNames(iChan);
imDir = MD.getChannelPaths(iChan);

im = imread([imDir{1} filesep imNames{1}{iFrame}]);

prot = MD.processes_{iProt}.loadChannelOutput;



%% ----- Figure Making ------ %%

imshow(im,[]);
figHan = gcf;
hold on

frameCols = colorMap(nFrames);

for j = 1:nFrames
   plot(prot.smoothedEdge{j}(:,1),prot.smoothedEdge{j}(:,2),'Color',frameCols(j,:));
end


