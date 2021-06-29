function [] = showTracksColorCoded(MD,varargin)
%function [] = showTracksColorCoded(MD) show movie overlaid with tracks
%color-coded with a property of choice
% input
%       MD                      MovieData file that has run FAPackage
%       Property                Names of a property of the tracks. E.g.:
%                               earlyAmpSlope
%       PropRange               The range of property value for min and max
%                               of the color-code
%       Colormap                The color map of choice (default: jet)
% output
%       a figure showing the movie and the overlay 
% Example: 
% showTracksColorCoded(MD,'Property','earlyAmpSlope','PropRange', [-300 1000]) 
% Sangyoon Han, May 31 2021

%% Input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('MD', @(x) isa(x,'MovieData'));
ip.addParameter('Property','earlyAmpSlope', @ischar); %0 for all classes and pick manually
ip.addParameter('PropRange',[-1 1], @isnumeric);
ip.addParameter('Colormap',jet, @(x) isnumeric(x) && size(x,2)==3); 
ip.parse(MD,varargin{:});

Property=ip.Results.Property;
PropRange=ip.Results.PropRange;
Colormap=ip.Results.Colormap;

tInterval = MD.timeInterval_;

%% Load FA package
faPackage=MD.getPackage(MD.getPackageIndex('FocalAdhesionPackage'));
% Load classification process
classProc = faPackage.getProcess(8);
iChan = find(classProc.checkChannelOutput);

%% persistent set up for large memory-requiring variable
persistent imgStack tracksNA curChanPath
%% Load tracksNA
finalProc = faPackage.getProcess(11);

%% Load imgStack, forceStack and anyother stack if it exists.
if isempty(curChanPath) || ~strcmp(curChanPath, MD.channels_(1).channelPath_) ...
    || isempty(tracksNA)
    tracksNA=finalProc.loadChannelOutput(iChan,'output','tracksNA');
end

if isempty(imgStack) || ~strcmp(curChanPath, MD.channels_(1).channelPath_)
    [imgStack] = getAnyStacks(MD);
    curChanPath = MD.channels_(1).channelPath_;
end
%% Launch pickAdhesion window with labeled adhesions with a right color and
% unlabed ones with white color. Get the right classes per newly selected
% adhesions
pickAdhesionTracksInteractive(tracksNA, imgStack,...
    'movieData',MD,'drawingOnly',true,'PropRange',PropRange,'Property',Property,'Colormap',Colormap);




