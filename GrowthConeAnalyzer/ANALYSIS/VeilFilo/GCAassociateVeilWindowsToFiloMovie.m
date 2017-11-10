function [ output_args ] = GCAassociateVeilWindowsToFiloMovie(movieData,varargin)
%DESCRIPTION: movieWrapper for marking 
% 
%  Input:
%
%   movieData - The MovieData object describing the movie, as created using
%   setupMovieDataGUI.m
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below:
%
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('OutputDirectory' -> character string) Optional. A character
%       string specifying the directory to save the backboneInfo structure to.
%       If not input, the backboneInfo will be saved in the same directory
%       as the movieData, in a sub-directory called "neurite_orientation_estimation_fixes"
%
%       ('ChannelIndex'-> Positive integer scalar or vector) The integer
%       indices of the channel(s) on which to perform the neurite orientation est.
%       This index corresponds to the channel's location in the array
%       movieData.channels_. If not input, all channels will be analyzed
%
%       ('ProcessIndex' -> Positive integer scalar or vector)
%       This parameter specifies the output process to use for performing the
%       estimation
%       This allows a previously processed image (ie shade corrected or
%       background subtracted to potentially be input). If not input, the
%       backbone information will be calculated from the channels (ie raw
%       images)
%
%%
% for now check movieData separately.
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end
%%Input check
ip = inputParser;

ip.CaseSensitive = false;

ip.addParameter('ChannelIndex',1);  

%ip.addParameter('TSOverlays',true);
defaultInDir = [movieData.outputDirectory_ filesep 'SegmentationPackage' filesep 'StepsToReconstruct'... 
    filesep 'VII_filopodiaBranch_fits_new' ];

ip.addParameter('InputDirectory',defaultInDir); % 

defaultOutDir = [movieData.outputDirectory_ filesep 'SegmentationPackage' filesep 'StepsToReconstruct' ... 
    filesep 'VIII_filopodiaBranch_WithVeil'];

ip.addParameter('OutputDirectory', defaultOutDir);

ip.addParameter('windMethod','ConstantNumber');

% ip.addParameter('OverlayColorByVel',true ); % this option will find the first and last window
% ip.addParameter('OverlayWindTracker',true); 
% ip.addParameter('OverlayOutlierFilter',true);  % outlier  
% ip.addParameter('OverlaySignalDetectType',{'P','R'}); % will only plot protrusion or retraction or quiescent events 

ip.addParameter('ReInit',61); 
% number to be defined in the subregion throughout the whole movie to
% define which windows are included. This allows for temporal gaps in the
% windows (ie where a window disappears and reappears) to be included in
% the subRoi.
ip.parse(varargin{:});

%% Init:
% nFrames = movieData.nFrames_;
nChan = numel(ip.Results.ChannelIndex);
imSize = movieData.imSize_;
ny = imSize(1); 
nx = imSize(2); 
 
for iCh = 1:nChan
  % INSERT CHECKS TO MAKE SURE THE FILO RECONSTRUCT WAS RUN  
  % load analInfo this contains the filoInfo for each frame    
  %load([p.OutputDirectory filesep ... 
   %   'Filopodia_Reconstruct_Channel_' num2str(iCh)' filesep 'analInfoTestSave.mat']); 
 % INSERT CHECKS TO MAKE SURE THE WINDOWING WAS RUN- Likely a nice way to
 % do this with  movieData to get all files associated with the windowing
 % step- for now
  load([ ip.Results.InputDirectory filesep 'Channel_' num2str(iCh) filesep 'filoBranch.mat']); 
  
    display(['Assigning Veil Windows To Filopodia for Channel' num2str(iCh)]); 
 
    if ~isdir(ip.Results.OutputDirectory)
        mkdir(ip.Results.OutputDirectory)
    end 
 
    %mkClrDir(saveDir)
    
    % get the list of image filenames
%     if p.ProcessIndex == 0
        imDir = movieData.channels_(iCh).channelPath_;
%     else
%         imDir = movieData.proceses_{p.ProcessIndex}.outfilePaths_;
%     end
    
    listOfImages = searchFiles('.tif',[],imDir,0);
    if isempty(listOfImages)
        error('No Images Found: Check Input Directory');  
    end
    
    
    % load the windowing directory : for now assume only one in the MD
    idxWind = cellfun(@(x) strcmpi(x.name_,'Windowing'),movieData.processes_);
    windDir  = movieData.processes_{idxWind}.outFilePaths_;
    windFiles = searchFiles('.mat',[],windDir,0,'all',1); 
    
%% %%%% Start Frame Loop %%%%
nFrames = length(filoBranch)-1; 
    for iFrame = 1:nFrames
       load(windFiles{iFrame})  ;    
       filoInfoC = filoBranch(iFrame).filoInfo; 
       filoInfo = GCAassociateVeilWindowsToFilo(filoInfoC,windows,0); 
       filoBranch(iFrame).filoInfo = filoInfo; % replace        
    end %frame loop 

% 
  save([ip.Results.OutputDirectory filesep 'filoBranchWindInfo.mat'],'filoBranch') % change name

end % channel loop 
end 

