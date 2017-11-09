function varargout = makeWindowClusterMovie(MD,varargin)
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
%       MD : MovieData object for the movie to make figure from.
%                   Optional. If not input, user will be asked to select
%                   a movieData file.
%
%       frameRange: Vector of frames to be displayed IMAGE from.
%                   
%       iChan : channel number 
%
%       winCluster : Cell array with the window clusters. Each cluster will be plotted with a different color 
%
%       colorMap : Handle to colormap to use for display, e.g. @jet or @hot
%                   Optional. Default is to use @jet.
%
%
%       framePath : This is the folder path in case you want to save the frames and make a movie.
%       
% Marco Vilela, 2013


%% ------------ Input ----------- %%

ip = inputParser;
ip.addRequired('MD',@(x) isa(x,'MovieData') );

%% ----------------- Init ---------- %%

winIdx  = MD.getProcessIndex('WindowingProcess');

%_________________________________________

ip.addParamValue('frameRange',1:MD.nFrames_);
ip.addParamValue('iChan',1,@isscalar);
ip.addParamValue('winCluster',{MD.processes_{winIdx}.nSliceMax_},@iscell);
ip.addParamValue('colorMap',@jet);
ip.addParamValue('framePath',[]);

ip.parse(MD,varargin{:});
frameRange = ip.Results.frameRange;
iChan      = ip.Results.iChan;
winCluster = ip.Results.winCluster;
colorMap   = ip.Results.colorMap;
framePath  = ip.Results.framePath;

%% ----- Figure Making ------ %%

nFramesSel = numel(frameRange);
nCluster   = numel(winCluster);
winColor   = num2cell(colorMap(nCluster),2);

fmt = ['%0' num2str(ceil(log10(nFramesSel))) 'd'];
ext = '.png';
zoom = 1;

ny = MD.imSize_(1);
nx = MD.imSize_(2);


% Configure figure

h = figure('Visible', 'on', 'Position', [50 50 nx ny]);

iptsetpref('ImshowBorder','tight');

set(h, 'InvertHardcopy', 'off');

set(h, 'PaperUnits', 'Points');

set(h, 'PaperSize', [nx ny]);

set(h, 'PaperPosition', [0 0 nx ny]); % very important

set(h,'DefaultLineLineSmoothing','on');

set(h,'DefaultPatchLineSmoothing','on');

axes('Position',[0 0 1 1]) 
hold on


for iFrame = frameRange
    
   iWindow = MD.processes_{winIdx}.loadChannelOutput(iFrame,iChan) ;
   imagesc(1:nx,1:ny,MD.channels_(iChan).loadImage(iFrame))
   axis([1 nx 1 ny])
   colormap(gray)
   hold on
   cellfun(@(x,y) plotWindows(iWindow(x),{y,'FaceAlpha',.2},0,'bandMin',1,'bandMax',1),winCluster,winColor','Unif',0);
   hold off
   
   if ~isempty(framePath)
       print(gcf, '-dpng', '-loose', ['-r' num2str(zoom*72)], [framesPath 'frame' num2str(iFrame, fmt) ext]);
   else
       pause(0.01)
   end
   
   cla(gca)
   
end

 
% Use these lines to make high res movies

%cmd = ['ffmpeg -loglevel quiet -y -r 15'  ' -i ' framesPath 'frame' fmt ext ' -b 50000k -bt 20000k ' framesPath 'movie.mp4 > /dev/null 2>&1'];
%system(cmd);


if nargout > 0
    varargout{1} = figHan;
end