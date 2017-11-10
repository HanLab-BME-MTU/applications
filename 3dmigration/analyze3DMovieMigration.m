function movieArray = analyze3DMovieMigration(movieArray,paramIn)
%ANALYZE3DMOVIEMIGRATION runs all processing and analysis steps for 3d single-cell migration
% 
% movieArray = analyze3DMovieMigration(movieArray)
% movieArray = analyze3DMovieMigration(movieArray,paramIn)
% 
% This function runs several processing and analysis steps to analyse the
% migration and branching behavior of the cell(s) in the input movie(s).
% This includes segmentation, skeletonization, pruning, tracking, etc. This
% is the primary analysis function of the 3dmigration project. Each
% processing and analysis step will be logged as a seperate process in the
% movieData object and stored in a separate folder of the analysis
% directory.
% 
% 
% 
% 
% Input:
% 
%   movieArray - A single MovieData3D object, or an array of MovieData3D
%   objects, also called a Movie Array.
% 
%   paramIn - A structure containing parameters to use in analysis. The
%   parameters/options available for this function are:
%
%       (FieldName->value)
%       
%       (RunParallel->Integer Scalar) If this number is greater than 1,
%       and more than one movieData was input, then the movies will be run
%       in a distributed fashion with one movie per processor, where the
%       number of processors is equal to RunParallel. If enabled, the
%       BatchMode option is also automatically enabled. The maximum value
%       for local parallelization is RunParallel = 8.
%       Optional. Default is 1 (no parallelization).
%
%       (BatchMode->true/false) If true, all graphical output is
%       suppressed, including progress bars, figures, input dialogues etc.
%       Optional. Default is false.
%
%       Parameters may also be input for individual processing steps. These
%       parameters should be in sub-fields of a field of the structure
%       which is named after the process the options are for. The available
%       options for each step are described in the help section for each
%       process function. The same parameters will be used for all movies
%       processed if a movieArray is input.
%
%           The process names and functions used are:
% 
%           
% 
%
% 
% 
% Hunter Elliott
% 12/2010
%

%% --------------- Input ---------------- %%

if nargin < 1 || isempty(movieArray) || ~isa(movieArray(1),'MovieData3D')
    error('The first input must be either a MovieData3D object or an array of MovieData3D objects!')
end

if nargin < 2 || isempty(paramIn)
    paramIn = [];
elseif ~isstruct(paramIn)
    error('The second input, paramIn, must be a structure!')
end



%% --------------- Init --------------- %%

nMov = numel(movieArray);


%List of processes to run
procList = {'SegmentationProcess3d',...
            'SkeletonizationProcess',...
            'SurfaceCurvatureProcess',...
            'BranchPruningProcess',...
            'BranchTrackingProcess',...
            
    



%% ------------ Processing ------------- %%
















