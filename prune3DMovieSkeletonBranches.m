function movieData = prune3DMovieSkeletonBranches(movieData,paramsIn)
%PRUNE3DMOVIESKELETONBRANCHES prunes some branches in the skeletons of the input movie using pruneSkeletonGraph.m
% 
% movieData = prune3DMovieSkeletonBranches(movieData3D)
% movieData = prune3DMovieSkeletonBranches(movieData3D,paramsIn)
% 
% This function calls pruneSkeletonGraph.m on the skeleton graphs of every
% frame of the input 3D movie and writes the pruned branches to disk in the
% specified output directory. The skeleton graphs for each frame should
% already have been created using skeletonize3DMovieMasks.m with the
% GetGraph option enabled. Additionally, this function uses the mask
% geometry for pruning, so the mask geometry must already have been
% analyzed using analyze3DMovieMaskGeometry.m
%
% Note: The original and pruned skeleton graphs are stored in image
% coordinates, which does not take into account the assymmetric pixel
% sizes. However, this function corrects for this fact prior to pruning so
% that the geometric pruning criteria are valid in all dimensions.
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
%       **NOTE** All of the options available in pruneSkeletonGraph.m can
%       also be set. They should be passed with a field name equal to the
%       option name in pruneSkeletonGraph.m
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
%% ------------------------- Parameters ------------------------%%

skName = 'pruned_skeleton_graph_frame_';%String for naming pruned skeleton files.

%% -------------------------- Input ---------------------------- %%

if nargin < 1 || isempty(movieData) || ~isa(movieData,'MovieData3D')
    error('The first input must be a valid MovieData3D object!');
end

if nargin < 2
    paramsIn = [];
end

%Get the indices of any previous pruning processes from this
%function
iProc = movieData.getProcessIndex('SkeletonPruningProcess',1,0);

%If the process doesn't exist, create it with default settings.
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(SkeletonPruningProcess(movieData));                                                                                                 
end

%Parse input, store in parameter structure. We also need to get the
%additional options for passing to pruneSkeletonGraph.m
[p,pruneParam] = parseProcessParams(movieData.processes_{iProc},paramsIn,true);

%% ------------------------ Init ------------------------------ %%


%Make sure the movie has been skeletonized, and that skeleton graphs were
%created.
iSkelProc = movieData.getProcessIndex('SkeletonizationProcess',1,~p.BatchMode);
if isempty(iSkelProc) 
    error('No skeletonization process found!')
else
    if isempty(p.ChannelIndex)    
        p.ChannelIndex = find(movieData.processes_{iSkelProc}.checkChannelSkeletonGraphs,1,'first');%There shouldn't be more than one skeleton channel, but just in case only find one
    end
    p.SkelProcessIndex = iSkelProc;    
    if isempty(p.ChannelIndex)
        error('No channels have valid skeleton graphs! Please create skeleton graphs first!')
    end        
end

%Make sure the mask geometry analysis has been run.
iMgProc = movieData.getProcessIndex('MaskGeometry3DProcess',1,~p.BatchMode);
if isempty(iMgProc) || ~movieData.processes_{iMgProc}.checkChannelOutput(p.ChannelIndex)
    error('No valid mask geometry analysis found! Please run mask geometry analysis first!')    
elseif movieData.processes_{iMgProc}.funParams_.PhysicalUnits
    %Lazy way to avoid re-scaling all the properties here.
    error('Please run mask geometry analysis with the physical units option set to false!')
end
p.MaskGeoProcessIndex = iMgProc;

%Get the segmentation process and verify that the masks are still good.
iSegProc = movieData.getProcessIndex('SegmentationProcess3D',1,~p.BatchMode);
if isempty(iSegProc) || ~movieData.processes_{iSegProc}.checkChannelOutput(p.ChannelIndex)
    error('The input MovieData does not have valid masks for the skeletonized channel!');
end
p.SegProcessInidex = iSegProc;

%Get mask file names and directory
maskDir = movieData.processes_{iSegProc}.outFilePaths_{1,p.ChannelIndex};
maskNames = movieData.processes_{iSegProc}.getOutMaskFileNames(p.ChannelIndex);

nFrames = movieData.nFrames_;

%Format string for zero-padding file names
fString = ['%0' num2str(floor(log10(nFrames))+1) '.f'];

%Set up the output directory
outDir = p.OutputDirectory;
mkClrDir(outDir);

if ~p.BatchMode
    wtBar = waitbar(0,'Please wait, pruning skeleton branches...');
end        

%Get the xy and z pixel sizes for scaling
if ~isempty(movieData.pixelSize_) && ~isempty(movieData.zSpacing_)
    pixXY = movieData.pixelSize_;
    pixZ = movieData.zSpacing_;
    hasSizes = true;    
else
    %warn the user, and assume unit pixel sizes.
    warning('Migration3D:MissingVoxelDimensions',...
        'Pixel XY size and Z spacing not specified in input movieData! Branch pruning will assume symmetric voxels of size 1nm!');
    pixXY = 1;
    pixZ = 1;
    hasSizes = false;
end

scFact = pixZ/pixXY;%The factor by which the z axis needs to be scaled to give symmetric voxels


%% ----------------------- Pruning -------------------------%%
%Go through each frame, prune the skeleton and save it to disk.

disp('Starting skeleton pruning...')


for iFrame = 1:nFrames

    
    %Load all the necessary files
    currMask = tif3Dread([maskDir filesep maskNames{1}{iFrame}]);
    currMG = movieData.processes_{iMgProc}.loadChannelOutput(p.ChannelIndex,iFrame);
    currSkGr = movieData.processes_{iSkelProc}.loadSkeletonGraph(p.ChannelIndex,iFrame);
    
    
    %----- Scaling ---- %
    %If voxel dimensions are available, scale the various inputs to take
    %this into account.
    if hasSizes && pixXY ~= pixZ
        
        %Scale the mask based on the pixel aspect ratio so the geometric
        %criteria used for pruning are correct.
        currMask = make3DImageVoxelsSymmetric(currMask,pixXY,pixZ);


        %Scale the skeleton branch coordinates to match the new symmetric voxels
        currSkGr.vertices(:,3) = currSkGr.vertices(:,3) .* scFact;
        hasPath = ~cellfun(@isempty,currSkGr.edgePaths);
        currSkGr.edgePaths(hasPath) = cellfun(@(x)([x(:,1:2) (x(:,3) .* scFact)]),...
                                       currSkGr.edgePaths(hasPath),'UniformOutput',false);


        %NOTE: The mask geometry properties already take the pixel aspect
        %ratio into account, so we don't need to scale them here.

    end

    % ----- Pruning ----- %
    
    %Call the pruning function on the properly scaled data, passing any
    %additional user-input pruning parameters.
    
    %We only prune if there is more than one edge! This keeps one bad frame
    %(with only a single edge) from crashing the entire process.
    if size(currSkGr.edges,1)> 1
        [vertices,edges,edgePaths,edgeLabels] = pruneSkeletonGraph(currSkGr.vertices,...
                                                        currSkGr.edges,...
                                                        currSkGr.edgePaths,...
                                                        currMask,...
                                                        currMG,...
                                                        pruneParam{:}); %#ok<NASGU,ASGLU>
    else
        vertices = currSkGr.vertices;
        edges = currSkGr.edges;
        edgePaths = currSkGr.edgePaths;
        edgeLabels = 2;
    end
    
    % ----- Output ----- %

    %Write the pruned skeleton to disk        
    numStr = num2str(iFrame,fString); %zero-pad the frame number    
    save([outDir filesep skName numStr '.mat'],'vertices','edges','edgePaths','edgeLabels')
    
    if ~p.BatchMode
        %Update the waitbar occasionally to minimize slowdown
        waitbar(iFrame/nFrames,wtBar)
    end
        
    
end



%% ---------------------- Finalization ---------------------------%%


%Store the input/output directories in the movieData
movieData.processes_{iProc}.setOutFilePath(p.ChannelIndex,outDir);
p.PruneParam = pruneParam;%Store the additional, function-determined parameters in the process array.
movieData.processes_{iProc}.setPara(p);
movieData.processes_{iProc}.setDateTime;
movieData.save;

if ~p.BatchMode && ishandle(wtBar)
    close(wtBar)
end

disp('Finished skeleton pruning!')


% TEMP - SAVE NEW PARAMS BECAUSE THEY HAVE BEEN ALTERED