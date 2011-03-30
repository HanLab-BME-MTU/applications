function movieData = prune3DMovieSkeletonBranches(movieData,paramIn)
%PRUNE3DMOVIESKELETONBRANCHES prunes some branches in the skeletons of the input movie using pruneSkeletonGraph.m
% 
% movieData = prune3DMovieSkeletonBranches(movieData3D)
% movieData = prune3DMovieSkeletonBranches(movieData3D,paramIn)
% 
% This function calls pruneSkeletonGraph.m on the skeleton graphs of every
% frame of the input 3D movie and writes the pruned branches to disk in the
% specified output directory. The skeleton graphs for each frame should
% already have been created using skeletonize3DMovieMasks.m with the
% GetGraph option enabled.
% 
% 
% Input: 
% 
%   movieData3D - The MovieData3D object describing the movie to prune
%   skeletons from. The movie must have already had skeleton graphs
%   created using skeletonize3DMovieMasks.m
% 
%   paramsIn - A structure containing the optional parameters to use. The
%   possible parameter field names and values are:
% 
%       (FieldName->fieldValue)
%
%       ('OutputDirectory' -> character string) Optional. A character
%       string specifying the directory to save the output to.
%       Optional. Default is a sub-directory of the movie's outputDirectory
%       called "pruned_branches"
%
%       (BatchMode->true/false) If true, all graphical output is
%       suppressed, such as progress bars, figures etc..
%
%   
% 
% Output: 
% 
%   movieData3D - The modified MovieData3D object with the parameters and
%   output locations logged in the processes_ array.
% 
%   Additionally, the pruned skeleton graphs will be written to a
%   sub-directory of the output directory.
% 
% 
% Hunter Elliott
% 3/2011
% 

%% -------------------------- Input ---------------------------- %%

if nargin < 1 || isempty(movieData) || ~isa(movieData,'MovieData3D')
    error('The first input must be a valid MovieData3D object!');
end

if nargin < 2
    paramIn = [];
end






%% ------------------------ Init ------------------------------ %%


%NEED TO SCALE EVERYTHING BY PIXEL ASPECT RATIO SO GEOMETRIC CRITERIA ARE
%CORRECT!!!!!!



