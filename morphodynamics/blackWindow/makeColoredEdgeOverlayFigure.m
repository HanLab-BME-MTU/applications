function varargout = makeColoredEdgeOverlayFigure(MD,iFrame,iEdgeFrames,iChan,colorMap)
%MAKECOLOREDEDGEOVERLAYFIGURE makes the figure with color-coded edge timelapse overlaid on a single image
% 
% makeColoredEdgeOverlayFigure;
% figHa = makeColoredEdgeOverlayFigure(movieData,imageFrameNumber,channelNumber,@colorMap);
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
% imageFrameNumber: Frame number to display IMAGE from.
%                   Optional. If not input, user will be asked.
%
% edgeFrameNumbers: Frame number(s) to overlay EDGES from.
%                   Optional. IF not input, user will be asked.
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
    iFrame = str2double(inputdlg('Enter frame number for the image to overlay edges on:','Image Frame Number'));
end

if isempty(iFrame)
    return
end

nFrames = MD.nFrames_;
if iFrame > nFrames || iFrame <= 0 || isnan(iFrame)
    error(['You must select a valid frame number to confinue!! Frame number must be between 1 and ' num2str(MD.nFrames_) '!!'])
end

if nargin < 3 || isempty(iEdgeFrames)
    [iEdgeFrames,OK] = listdlg('ListString',arrayfun(@num2str,1:nFrames,'Unif',false), ...
        'PromptString','Select edge frame(s) to overlay:','ListSize',[300 400]);   

    if ~OK
        return
    end
end



nChan = numel(MD.channels_);
if nChan == 1
    iChan = 1;
    OK = true;
elseif nargin < 4 || isempty(iChan)
    [iChan OK] = listdlg('ListString',arrayfun(@num2str,1:nChan,'Unif',false),'SelectionMode','single',...
        'PromptString','Select a channel to display image from:','ListSize',[300 400]);
end
if ~OK
    return
end

if nargin < 5 || isempty(colorMap)
    colorMap = @jet;
end


%% ----------------- Init ---------- %%

imNames = MD.getImageFileNames(iChan);
imDir = MD.getChannelPaths(iChan);
im = imread([imDir{1} filesep imNames{1}{iFrame}]);
prot = MD.processes_{iProt}.loadChannelOutput;


%% ----- Figure Making ------ %%

figHan = figure;
imshow(im,[]);
hold on
nFramesSel = numel(iEdgeFrames);
frameCols = colorMap(nFramesSel);
for j = 1:nFramesSel
   plot(prot.smoothedEdge{iEdgeFrames(j)}(:,1),prot.smoothedEdge{iEdgeFrames(j)}(:,2),'Color',frameCols(j,:));
end


if nargout > 0
    varargout{1} = figHan;
end

