function [ output_args ] = GCAVisualsProtrusionMapMovie(movieData,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% INPUT
%
%   movieData - The MovieData object describing the movie, as created using
%   movieSelectorGUI.m
%
%   Parameter Structure Field Names:
%
% Generic Fields: (Input/Output Fields Needed for Wrapper)
%       ('OutputDirectory' -> Optional. A character
%       string specifying the directory to save the Visualization Output
% for now check movieData separately.
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end

% load the protrusionSamples 

%% Check Input 
%cDir = [MD.outputDirectory_ filesep 'protrusion_samples'];
ip = inputParser;

ip.CaseSensitive = false;

ip.addParameter('OutputDirectory', pwd);
ip.addParameter('CLims',[-100,100]); % currently in nm/sec
ip.addParameter('smooth',false); 
ip.parse(varargin{:});

%%
idxSampProc =  cellfun (@(x) strcmpi(x.name_,'Protrusion Sampling'),movieData.processes_);
cFile = movieData.processes_{idxSampProc}.outFilePaths_; 

load(cFile{1}); 
 mapValues =  protSamples.avgNormal;
 mapValues = mapValues*movieData.pixelSize_/movieData.timeInterval_;
cmap = brewermap(128,'RdBu');
cmap = flip(cmap,1);
if ip.Results.smooth
    mapValues= smoothActivityMap(mapValues,'upSample',1,'SmoothParam',0.75);
end
GCAVisualsProtrusionMap(mapValues,'CLim',ip.Results.CLims,'colorMap',cmap,'frameInSec',movieData.timeInterval_); 

saveas(gcf,[ip.Results.OutputDirectory filesep 'protrusionMap.fig']); 
saveas(gcf,[ip.Results.OutputDirectory filesep 'protrusionMap.eps']); 
helperScreen2png([ip.Results.OutputDirectory filesep 'protrusionMap.png']); 

close gcf
end

