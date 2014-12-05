function [ output_args ] = GCAassociateVeilWindowsToFiloMovie(movieData)
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
%% CHECK Parameters 
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end

if nargin < 2
    % Generic
    paramsIn.OutputDirectory = [movieData.outputDirectory_ filesep 'filopodia_fits'];
    paramsIn.ChannelIndex = 1;
    paramsIn.ProcessIndex = 0; % use raw images
end 
p = paramsIn;

%% Init:
nFrames = movieData.nFrames_;
nChan = numel(p.ChannelIndex);
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
  load([p.OutputDirectory filesep ... 
      'Filopodia_Fits_Channel_' num2str(iCh) filesep 'analInfoTestSave.mat']); 
  
    display(['Assigning Veil Windows To Filopodia for Channel' num2str(iCh)]); 
 
   
    % make final output dir where Assignment Info will be saved 
    saveDir =  [p.OutputDirectory  filesep ...
        'Filopodia_Fits_Channel_'  num2str(iCh)]; % I think just save back 
    % into the same directory for now. 
    %mkClrDir(saveDir)
    
    % get the list of image filenames
    if p.ProcessIndex == 0
        imDir = movieData.channels_(iCh).channelPath_;
    else
        imDir = movieData.proceses_{p.ProcessIndex}.outfilePaths_;
    end
    
    listOfImages = searchFiles('.tif',[],imDir,0);
    if isempty(listOfImages)
        error('No Images Found: Check Input Directory');  
    end
%% %%%% Start Frame Loop %%%%
    for iFrame = 1:nFrames
       load([movieData.outputDirectory_ filesep 'windows' filesep  'windows_frame__frame_' num2str(iFrame,'%03d') '.mat'])  ;    
       filoInfoC = analInfo(iFrame).filoInfo; 
       filoInfo = GCAassociateVeilWindowsToFilo(filoInfoC,windows,0); 
       analInfo(iFrame).filoInfo = filoInfo; % replace 
         
        
    end %frame loop 

% 
  save([saveDir filesep 'analInfoWithWindInfo.mat'],'analInfo') % change name




end % channel loop 
end 

