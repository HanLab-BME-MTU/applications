function movieData = skeletonize3dMovie(movieData,paramsIn)
%SKELETONIZEMOVIE creates and saves skeletons for the masks in every frame of the input movie
% 
% movieData3D = skeletonize3dMovie(movieData3D)
% movieData3D = skeletonize3dMovie(movieData3D,paramsIn)
% 
% This function skeletonizes the masks for every frame of the 3D movie
% described by the input MovieData3D object. The skeleton is created by
% thinning. The skeleton is also broken into labelled edges and vertices,
% with each separate edge and vertex labelled with a different integer
% value in the output matrix, similar to how the bwlabel.m function labels
% objects in masks. Additionally, the ordered X,Y,Z coordinates of each
% edge can be calculated and output. This additional ourput will slow the
% processing.
% Additionally, the skeletons for each frame may be converted to a graph
% structure, using skel2graph.m, and these outputs saved. This will
% significantly increase processing time.
% 
%   ***WARNING*** Be careful when specifying the output directory, as any
%   existing files in this directory will be deleted to make room for the
%   skeletons.
%
% NOTE: This function uses the thinning-based algorith implemented in
% skeleton3d.m
% 
% Input:
% 
%   movieData3D - A MovieData3D object describing the movie to skeletonize.
%   This movie must already have been segmented, and have a valid
%   SegmentationProcess3D in its process array.
% 
%   paramsIn - A structure containing the optional parameters to use. The
%   possible parameter field names and values are:
% 
%       (FieldName->fieldValue)
%
%       ('GetGraph'->true/false) If true, the graph structure of the
%       skeleton will be calculated using skel2graph.m and saved to disk.
%       Enabling this option will significantly increase the time required
%       for this function to run.
%       Optional. Default is false.
%
%       ('ClearBoundary'->true/false) If true, pixels on the image boundary
%       will be removed prior to skeletonization. These pixels will not be
%       eroded during the skeletonization by the method this function uses.
%       Optional. Default is true.
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
%   NOTE: The skeleton files will be bit-packed with the default matlab
%   compression system, because they are generally very sparse. This may
%   preclude opening the files with ImageJ or other software.
%
% Hunter Elliott
% 12/2010
%

%%  ------------------ Input ------------------ %%

if nargin < 1 || isempty(movieData) || ~isa(movieData,'MovieData3D')
    error('The first input must be a valid MovieData3D object!')
end

if nargin < 2 || isempty(paramsIn)
    paramsIn = [];
end

%Get the indices of any previous skeletonization processes from this
%function
iProc = movieData.getProcessIndex('SkeletonizationProcess',1,0);

%If the process doesn't exist, create it with default settings.
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(SkeletonizationProcess(movieData));                                                                                                 
end

%Parse input, store in parameter structure
p = parseProcessParams(movieData.processes_{iProc},paramsIn);

% ---- Parameter Checks ---- %

if numel(p.ChannelIndex) ~= 1
    error('You must select exactly 1 channel to use masks from for skeletonization!')
end


%% ------------------ Init -------------------- %%


%Make sure the movie has been segmented and has masks
iSegProc = movieData.getProcessIndex('SegmentationProcess3D',1,~p.BatchMode);

if isempty(iSegProc) || ~movieData.processes_{iSegProc}.checkChannelOutput(p.ChannelIndex)
    error('The input MovieData does not have valid masks for the selected channel!');
end

%Get mask file names and directory
maskDir = movieData.processes_{iSegProc}.outMaskPaths_{p.ChannelIndex};
maskNames = movieData.processes_{iSegProc}.getOutMaskFileNames(p.ChannelIndex);

%Set up output directories
mkClrDir(p.OutputDirectory);
skelDir = p.OutputDirectory;

%And separate directories for graphs and skeletons if requested
if p.GetGraph
    graphDir = [p.OutputDirectory filesep 'skeleton_graphs'];
    mkClrDir(p.OutputDirectory);
end

nFrames = movieData.nFrames_;
nSlices = movieData.nSlices_;

if ~p.BatchMode
    wtBar = waitbar(0,'Please wait, skeletonizing masks ...');
end        

%Format string for zero-padding file names
fString = ['%0' num2str(floor(log10(nFrames))+1) '.f'];



%% --------------- Skeletonization -------------------- %%

disp('Starting skeletonization...')
disp(['Using masks from ' maskDir ', results will be saved to ' skelDir ])

for iFrame = 1:nFrames
    
   
    %Load the current mask
    currMask = tif3Dread([maskDir filesep maskNames{1}{iFrame}]);
        
    
    %Remove boundary pixels from skeleton if requested.
    if p.ClearBoundary
        currMask(:,:,[1 end]) = false;
        currMask(:,[1 end],:) = false;
        currMask([1 end],:,:) = false;
    end
    
    %Skeletonize it
    currSkel = skeleton3D(currMask);
        
    %We want to compress the masks, so don't use stackWrite.m
    for i = 1:nSlices
        %Append each z-slice to the tiff
        imwrite(currSkel(:,:,i),[skelDir filesep ...
        'skeleton_' maskNames{1}{iFrame}(1:end-3) 'tif'],'tif','WriteMode','append')
    end
    
    if p.GetGraph
        %This skeletonization method produces 26-connected skeletons so we
        %need to use this connectivity for graph-ifying the skeleton.
        [vertices,edges,edgePaths] = skel2graph(currSkel,26);      %#ok<NASGU,ASGLU>   
          
        %Save graph structure to file    
        numStr = num2str(iFrame,fString); %zero-pad the frame number
        save([graphDir filesep 'skeleton_graph_frame_' numStr '.mat'],...
                'vertices','edges','edgePaths');
    end
  
    if ~p.BatchMode && mod(iFrame,5)
        %Update the waitbar occasionally to minimize slowdown
        waitbar(iFrame/nFrames,wtBar)
    end
    
end

%% ---------- Finalization ---------- %%


if ~p.BatchMode && ishandle(wtBar)
    close(wtBar)
end

movieData.processes_{iProc}.setDateTime;
movieData.saveMovieData;

disp('Finished Skeletonizing!')


