function movieData = skeletonize3dMovie(movieData,paramsIn)
%SKELETONIZEMOVIE creates and saves skeletons for the masks in every frame of the input movie
% 
% movieData3D = skeletonize3dMovie(movieData3D)
% movieData3D = skeletonize3dMovie(movieData3D,paramIn)
% 
% This function skeletonizes the masks for every frame of the 3D movie
% described by the input MovieData3D object. The skeleton is created by
% thinning. The skeleton is also broken into labelled edges and vertices,
% with each separate edge and vertex labelled with a different integer
% value in the output matrix, similar to how the bwlabel.m function labels
% objects in masks. Additionally, the ordered X,Y,Z coordinates of each
% edge can be calculated and output. This additional ourput will slow the
% processing.
% 
% 
% 
% Input:
% 
%   movieData3D - A MovieData3D object describing the movie to skeletonize.
%   This movie must already have been segmented, and have a valid
%   SegmentationProcess3D in its process array.
% 
%   paramIn - A structure containing the optional parameters to use. The
%   possible parameter field names and values are:
% 
%       (FieldName->fieldValue)
%
%       ('OutputDirectory' -> character string) Optional. A character
%       string specifying the directory to save the skeletons to.
%       Optional. Default is a sub-directory of the movie's outputDirectory
%       called "skeletonization"
%
%       ('ChannelIndex' -> positive integer) Integer indices of the
%       channel to use masks from for skeletonization.
%       Optional. Default is channel 1.
%
%       (BatchMode->true/false) If true, all graphical output is
%       suppressed, such as progress bars, figures etc..
%
%       
% Output:
%
%   movieData3D - The modified MovieData3D object with the skeletonization
%   logged in the processes array.
%
%   Additionally, the skeletons will be written to the selected
%   outputDirectory as multi-page .tif image files, one per frame of the
%   movie. If ordered edge path coordinates are returned, these will be
%   written to a sub-directory of the directory containing the skeleton
%   image files.
%
% Hunter Elliott
% 12/2010
%

%%  ------------------ Input ------------------ %%

if nargin < 1 || isempty(movieData) || ~isa(movieData,'MovieData3D')
    error('The first input must be a valid MovieData3D object!')
end

if nargin < 2 || isempty(paramIn)
    paramsIn = [];
end

%Get the indices of any previous skeletonization processes from this
%function
iProc = movieData.getProcessIndex('SkeletonizationProcess',1,0);

%If the process doesn't exist, create it with default settings.
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(SkeletonizationProcess(movieData,movieData.outputDirectory_));                                                                                                 
end

%Parse input, store in parameter structure
p = parseProcessParams(movieData.processes_{iProc},paramsIn);



%% ------------------ Init -------------------- %%


%Make sure the movie has been segmented and has masks
iSegProc = movieData.getProcessIndex('SegmentationProcess3D',1,~p.BatchMode);

if isempty(iSegProc) || ~movieData.processes_{iSegProc}.checkChannelOutput(p.ChannelIndex)
    error('The input MovieData does not have valid masks for the selected channel!');
end






