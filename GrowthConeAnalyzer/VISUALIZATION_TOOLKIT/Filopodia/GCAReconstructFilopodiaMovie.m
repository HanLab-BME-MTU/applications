function [ output_args ] = GCAReconstructFilopodiaMovie(movieData,paramsIn )
%
% GCAreconstructFilopodiaMovie.m :  movieWrapper function for
% GCAReconstructFilopodia.m 
% STEP V in the GCA package
%% INPUT:
% 
%   movieData - The MovieData object describing the movie, as created using
%   movieSelectorGUI.m 
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below:
% 
%   
%   Parameter Structure Field Names:
%   
% Generic Fields: (Input/Output Fields Needed for Wrapper) 
%       ('OutputDirectory' -> character string) Optional. A character
%       string specifying the directory to save the backboneInfo structure to.
%       If not input, the backboneInfo will be saved in the same directory 
%       as the movieData, in a sub-directory called
%       "neurite_orientation_estimations"  
%       % PERSONAL NOTE : You might need
%       to truncate these names for the windows paths. % 
%
%       ('ChannelIndex'-> Positive integer scalar or vector) The integer
%       indices of the channel(s) on which to perform the neurite orientation 
%       estimation. 
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
%  GCAReconstructFilopodia specific input
%
%       ('protrusion -> output of  
% 
%       
%  
% OUTPUT: (see main function GCAReconstructFilopodia- for details): 
%       adds fields to analInfo Nx1 structure array where N = the number of frames 
%       addedFields 



end

