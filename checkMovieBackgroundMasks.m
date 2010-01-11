function status = checkMovieBackgroundMasks(movieData,iChannels)
% 
% status = checkMovieBackgroundMasks(movieData,iChannels)
% 
% This function returns true if the movie has masks background masks for
% the selected channels and false otherwise.
%
% Input: 
%
%   movieData - The structure describing the movie created with
%   setupMovieData.m
%
%   iChannels - Optional. The integer indices of the channels to check
%   masks for. If not input, only the status of the masks as indicated in
%   the movieData structure will be checked (this is much faster, but
%   riskier).
%
%
% Output:
% 
%   status - Logical. True if masks have been successfully created, false
%   otherwise.
% 
%
% Hunter Elliott, 11/2009
%


if nargin < 2
    iChannels = [];
end
    

status = false;

if isfield(movieData,'backgroundMasks') && ...    
    isfield(movieData.backgroundMasks,'directory') && ...
    exist(movieData.backgroundMasks.directory,'dir') ...
    && isfield(movieData.backgroundMasks,'channelDirectory')...
    && movieData.backgroundMasks.status == 1
    
    
    if isempty(iChannels)
        status = true;
        return
    elseif any(iChannels > length(movieData.backgroundMasks.channelDirectory)) ...
            || any(iChannels < 1)
        status = false;
        return
    else
        
        %If channels were explicitly requested, check the channel
        %directories and masks
                
        nCheck = length(iChannels);
        status = false(1,nCheck);
        
        for i = 1:nCheck
            
            status(i) = ~isempty(movieData.backgroundMasks.channelDirectory{iChannels(i)}) && ...
                exist([movieData.backgroundMasks.directory filesep ...                
               movieData.backgroundMasks.channelDirectory{iChannels(i)}],'dir') && ...
                length(dir([movieData.backgroundMasks.directory filesep ...
                movieData.backgroundMasks.channelDirectory{iChannels(i)} filesep '*.tif'])) == movieData.nImages(iChannels(i));
    
    
        end
   end
end
    
status = all(status);
