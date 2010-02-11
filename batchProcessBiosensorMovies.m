function movieArray = batchProcessBiosensorMovies(movieArray,actChanName,volChanName,varargin)
%BATCHPROCESSBIOSENSORMOVIES processes a series of biosensor movies
%
% movieArray = batchProcessBiosensorMovies(movieArray,actChanName,volChanName)
% 
% movieArray = batchProcessBiosensorMovies(movieArray,actChanName,volChanName,'OptionName',optionvalue...)
% 
%
% This function produces activity ratio images for a set of input
% ratiometric biosensor movies, by calling the function
% processBiosensorMovie.m on each. All movies will be processed using the
% same options/parameters.
% 
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
%   'OptionName',optionValue - A string with an option name followed by the
%   value for that option. If these option names are not recognized by this
%   function, the options will be passed to the processing function for
%   each movie. See processBiosensorMovie.m for possible options.
% 
%   Possible Option Names:
%
%       ('OptionName' -> possible values)
%
%       ('ActivityShadeChanName' -> character string / regexp)
%       This option specifies the name of the channel containing the
%       shade-correction images for the activity channel.
%       Optional. If not input, the software will look for channels whose
%       name contains the activity channel name and the word "shade" If it
%       cannot find a channel, the movie will be skipped.
% 
%       ('VolumeShadeChanName' -> character string / regexp)
%       This option specifies the name of the channel containing the
%       shade-correction images for the volume channel.
%       Optional. If not input, the software will look for channels whose
%       name contains the volume channel name and the word "shade" If it
%       cannot find a channel, the movie will be skipped.
%
% Hunter Elliott
% 1/2010
%

%% ------ Parameters ----- %%

shStr = 'shade'; %The string to look for in the shade correction channel names.


%% ----------- Input -------- %%

if nargin < 3
    error('Must input a movieArray, and the name of the activity and volume channels! Check input!')
end

nMov = length(movieArray); % Get number of movies to process

if ~iscell(movieArray) || nMov < 1
    error('First input, movieArray, must be a cell-array containing at least one movieData structure!')
end

if ~(ischar(actChanName) && ischar(volChanName))
    error('Inputs 2 and 3 must be character strings with the names of the activity and volume channels respectively!')
end

%Parse the input, seperating inputs for this function from inputs for
%processing functions
[commonOptions,actShadeName,volShadeName] = parseInputs(varargin);

%If a name wasn't specified, create a regexp to search for correctly named
%channels for the shade correction images.
if isempty(actShadeName) 
    actShadeName = ['(' actChanName '.*' shStr ')|(' shStr '.*' actChanName ')'];
end
if isempty(volShadeName) 
    volShadeName = ['(' volChanName '.*' shStr ')|(' shStr '.*' volChanName ')'];
end



%% --------- Init ---------- %%

%Check if the batchMode setting has been specified, and if not set it
if ~any(strcmp('BatchMode',commonOptions))    
   commonOptions = [commonOptions {'BatchMode',true}]; %Default is to suppress output from sub-functions         
end


%% ------------ Processing -------------- %%
%Go through each movie and call the processing function.

for iMov = 1:nMov
            
    
    try
                
    
        %Get the number indices for the specified volume and activity channels,
        %as this may vary from movie-to-movie.
        iAct = find(strcmp(actChanName,movieArray{iMov}.channelDirectory),1); %Activity channel index
        iVol = find(strcmp(volChanName,movieArray{iMov}.channelDirectory),1); %Volume channel index               
        
        %Get indices for shade-correction images
        iActShade = find(cellfun(@(x)(~isempty(regexpi(x,actShadeName))),movieArray{iMov}.channelDirectory)); %Activity channel shade-correction images
        iVolShade = find(cellfun(@(x)(~isempty(regexpi(x,volShadeName))),movieArray{iMov}.channelDirectory)); %Volume channel shade-correction images
        
        %Make sure only one shade-correction channel was found
        if length(iActShade) ~= 1 || length(iVolShade) ~= 1
            error('Could not unambigously find shade-correction image channels! Retry specifying shade image channels!')                        
        end
        
        %Combine the selected channels for this movie with the options used
        %on all movies
        currOptions = [commonOptions {'shadeCorrection_ShadeImageChannels',[iActShade iVolShade], ...
            'ActivityChannel', iAct,'VolumeChannel',iVol}];
        
        %Process the current movie
        movieArray{iMov} = processBiosensorMovie(movieArray{iMov},currOptions{:});
        
        %If there was an error during processing of this movie previously,
        %remove it to reflect successful completion
        if isfield(movieArray{iMov}.biosensorProcessing,'error')
            movieArray{iMov}.biosensorProcessing = ...
                rmfield(movieArray{iMov}.biosensorProcessing,'error');            
        end
    
    catch errMess

        movieArray{iMov}.biosensorProcessing.error = errMess;
        
        disp(['Error processing movie ' num2str(iMov) ' : ' errMess.message]);            
    
    end
    
    %Save the modified movieData structure
    updateMovieData(movieArray{iMov})
    
    
end

function [commonOptions,actShadeName,volShadeName] = parseInputs(argArray)

%---Defaults---%
actShadeName = [];
volShadeName = [];
commonOptions = [];

if isempty(argArray)
    return
end

nArg = length(argArray);

commonOptions = cell(1,nArg);

%Make sure there is an even number of arguments corresponding to
%optionName/value pairs
if mod(nArg,2) ~= 0
    error('Inputs must be as optionName / value pairs!')
end

for i = 1:2:nArg

    switch argArray{i}

        case 'ActivityShadeChanName'
            
            actShadeName = argArray{i+1};
            
        case 'VolumeShadeChanName'
            
            volShadeName = argArray{i+1};
            
    
        otherwise
            %If the option isn't recognized, store it in an array so it can
            %be passed to the processing function.
            commonOptions(i:i+1) = argArray(i:i+1);


    end

end

%Remove empty entries in commonOptions
commonOptions = commonOptions(cellfun(@(x)(~isempty(x)),commonOptions));


