function [ output_args ] = performTrackingMovie(movieData,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% params 


% for now check movieData separately.
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end
%%Input check
ip = inputParser;

ip.CaseSensitive = false;

% PARAMETERS

ip.addParameter('InputDirectory',[]); 
ip.addParameter('ChannelIndex',1); 
ip.addParameter('frames',[]); 

ip.addParameter('OutputDirectory',[]); 

% CostMatParams 
ip.addParameter('searchRadius',5); % search radius to be used for tracking  
ip.addParameter('predict',false);  % whether or not to predict the filopodia's position based 
% on the protrusion vector measurement. 

ip.parse(varargin{:});
p = ip.Results;

%%
if isempty(ip.Results.InputDirectory)
    inDir =  [movieData.outputDirectory_ filesep 'SegmentationPackage' ...
        filesep 'StepsToReconstruct' filesep...
        'VII_filopodiaBranch_fits' ];
else
    inDir = ip.Results.InputDirectory;
end

if isempty(ip.Results.OutputDirectory)
    
    outDir = [movieData.outputDirectory_ filesep 'tracking'];
else
    outDir = ip.Results.OutputDirectory;
end
if ~isdir(outDir);
    mkdir(outDir);
end


costMatParams.searchRadius = ip.Results.searchRadius;
costMatParams.predict = ip.Results.predict;

%% 20140915 for now just load analInfo automatically from reconstruct channel 1
% and images from donor channel under the current output set-up

% load from
load([inDir filesep 'Channel_' num2str(ip.Results.ChannelIndex) filesep 'filoBranch.mat']);

imDir = movieData.channels_(1).channelPath_;

performTracking(filoBranch,costMatParams,imDir, outDir)

end


