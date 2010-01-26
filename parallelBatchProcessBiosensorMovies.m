function movieArray = parallelBatchProcessBiosensorMovies(movieArray,...
                                                          actChanName,...
                                                          volChanName,...
                                                          nProcessors,...
                                                          varargin)
                                                      
% PARALLELBATCHPROCESSBIOSENSORMOVIES parallel processes a group of biosensor movies
% 
% movieArray = parallelBatchProcessBiosensorMovies(movieArray,...
%                                                  actChanName,...
%                                                  volChanName,...
%                                                  nProcessors,...
%                                                  'OptionName',optionValue
%                                                  )                                              
%                                                       
%                                                       
% This function distributes the processing of the movies in the input
% movieArray to  run on multiple processors using the distributed computing
% toolbox to call batchProcessBiosensorMovies.m on sub-sets of the input
% movies.                                                       
%
%
% Input:
%                                                       
%   movieArray - A cell array of movieData structures describing the movies
%   to analyze. The individual movieData structures are as created with
%   setupMovieData.m. An array can be created from existing movieDatas
%   using setupMovieArray.m
%
%   actChanName - A character string containing the name of the channel
%   which contains the activity images (numerator of ratio). This is also
%   the name of the sub-directory within the image directory which contains
%   these images.
%
%   volChanName - A character string containing the name of the channel
%   which contains the volume images (denominator of ratio). This is also
%   the name of the sub-directory within the image directory which contains
%   these images.
%
%   nProcessors - A positive integer <= 8 which specifies the number of
%   processors to distribute the movie processing among.
%   Optional. Default is 8 (the maximum for local processing). 
%
%   'OptionName',optionValue - A string with an option name followed by the
%   Value for that option. See batchProcessBiosensorMovies.m for option
%   descriptions and possible values (all options will be passed to this
%   function).                                                     
%                                                       
%                                                       
%                                                       
% Output:                                                      
%                                                       
%  movieArray - The cell array of movieData structures will all processing
%  steps and parameters logged.                                                     
%                                                       
%                                                       
%                                                       
% Hunter Elliott                                                      
% 1/2010   

%% --------------- Input -------------- %%

if nargin < 4 || isempty(nProcessors)
    nProcessors = 8;
end

if nargin < 3
    error('Must input a movieArray, and the name of the activity and volume channels! Check input!')
end

nMov = length(movieArray); % Get number of movies to process

if ~iscell(movieArray) || nMov < 2
    error('First input, movieArray, must be a cell-array containing at least two movieData structures!')
end

if ~(ischar(actChanName) && ischar(volChanName))
    error('Inputs 2 and 3 must be character strings with the names of the activity and volume channels respectively!')
end


%% ---------- Init ---------- %%

disp('Initializing parallel processing...')

%Get the size of the current pool, if any
pSize = matlabpool('size'); 
if pSize ~= nProcessors
   
    %If one is open, close it
    if pSize > 0        
        matlabpool('close')            
    end        
    
    %Initialize the workers for this job
    matlabpool('open',nProcessors)    
    
end

%Slice the options for distribution
optionsIn = [varargin {'BatchMode', true}];
optionsIn = repmat(optionsIn,nMov,1);

%% ------- Processing ------ %%

parfor iMov = 1:nMov

    disp(['Processing movie ' num2str(iMov)]);

    movieArray(iMov) = batchProcessBiosensorMovies(movieArray(iMov),actChanName,volChanName,optionsIn{iMov,:});

end


disp('Finishing parallel processing...')

matlabpool('close')










