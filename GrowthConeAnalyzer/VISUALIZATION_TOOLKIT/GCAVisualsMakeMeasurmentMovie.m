function [ output_args ] = GCAVisualsMakeMeasurmentMovie(MD,varargin)
% GCAVisualsMakeMeasurementMovie 
% 
% 
% This function will make the measurement movies for all measurements 
% Default is to ask the user which type of movie they would like to make 

% INPUT: MD movieData object
% 
%% Check input
ip = inputParser;

ip.CaseSensitive = false;

ip.addParameter('interactive',true,@(x) islogical(x));

% Note need to make more generic for channel wrap! 
defaultMeasDir = [MD.outputDirectory_ filesep 'GCAMeasurementExtraction' filesep 'Channel_1']; 
ip.addParameter('MeasurementDirectory',defaultMeasDir,@(x) ischar(x)); 

defaultInDir = [MD.outputDirectory_ filesep 'GrowthConeAnalyzer' filesep 'SegmentationPackage' ... 
    filesep 'StepsToReconstruct' filesep 'VII_filopodiaBranch_fits']; 

ip.addParameter('InputDirectory', defaultInDir,@(x) ischar(x)); 

ip.parse(varargin{:});
p = ip.Results;

%% go into each folder and look for a measurement 
% for now just search files - redesign so that the parameters in the
% future are more cleverly named

% might also include ylabel name and ylim for each parameter and
% read in each Descriptor directory to keep constant.

% search all descriptor parameters.
localParamFiles = searchFiles('param_',[],[ip.Results.MeasurementDirectory filesep 'Descriptor'],1);

paramNamesC = cellfun(@(x) strrep(x,'param_',''),localParamFiles(:,1),'uniformoutput',0);
paramNamesC = cellfun(@(x) strrep(x,'.mat',''),paramNamesC,'uniformoutput',0);

% if interactive
% ask the user for which measurement they would like to make a sanity movie

% most measurments are going to be simple it will just be plotting
% the filopodia of the filter set with the value saved
if ip.Results.interactive == true
    paramSelect  = listSelectGUI(paramNamesC,[],'move');
    selected  = paramNamesC(paramSelect);
else % make all movies
    selected = paramNamesC; 
end

%% Different File types might require different functions
for iSelect = 1:numel(selected)  
    cParam = selected{iSelect}; 
%    switch  strcmpi(
% if ~isempty(regexpi('filoDensityAlongBranch',cParam)); 
    
    
end 


end

