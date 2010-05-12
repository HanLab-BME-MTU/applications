function movieData = refineMovieMasks(movieData,varargin)
% REFINEMOVIEMASKS Performs post-processing to improve masks for an input movie.
% 
% movieData = refineMovieMasks(movieData)
% 
% movieData = refineMovieMasks(movieData,'OptionName',optionValue)
% 
% 
% This function performs several post-processing steps to refine the masks
% for the input movie, overwriting the existing masks. The available
% post-processing steps are listed below.
% 
% 
% Input:
% 
%   movieData - The structure describing the movie, as created using
%   setupMovieData.m
% 
%   Possible Options:
%       ('OptionName' -> possible values)
% 
%
%       ('ChannelIndex'-> Positive integer scalar or vector) 
%       The integer indices of the channel(s) to perform background subtraction
%       on. This index corresponds to the channel directories location in
%       the cell array movieData.channelDirectory. If not input, the user
%       will be asked to select from the available channels
%
%       ('MaskCleanUp' -> True/False)
%       If true, various operations will be performed to clean up the
%       mask. The operations are listed below, and will be performed in the
%       order they are listed.
%       Optional. Default is True.
%
%       These operations include:
%
%           Clean-up Methods/Parameters:
%
%           ('MinimumSize' -> Positive integer scalar)
%           Any objects in the mask whose area is below this size (in
%           pixels) will be removed via morphological area opening.
%           Optional. Default is 10 pixels. Set to Inf to keep all objects.
%        
%           ('ClosureRadius' -> non-negative integer scalar)
%           If this is parameter greater than zero, the mask will be closed
%           using a disk-shaped structuring element of this radius. This
%           has the effect of connecting previously un-connected components
%           in the mask if they are within 2x this distance of one another.
%           Optional. Default is 3 pixels.
%
%           ('ObjectNumber -> Positive integer scalar)
%           Only this number of the largest objects in the mask will be
%           kept. That is, if this number is 2, only the two largest
%           objects will be kept. This step is performed AFTER the edge
%           refinement (if enabled)
%           Optional. Default is 1. Set to Inf to keep all objects.
%
%           ('FillHoles -> True/False)
%           If true, any holes in any objects in the mask will be filled
%           in.
%           Optional. Default is true.
%
%       ('EdgeRefinement' -> True/False)
%       If true, edge detection will be used to refine the position of the
%       mask edge to the nearest detected edge in the image. This will be
%       done after any of the cleanup procedures listed above.
%       Optional. Default is False.
% 
%           Edge Refinement Parameters:
%
%           NOTE: For descriptions and default values, see refineMaskEdges.m           
%           
%           ('MaxEdgeAdjust' -> Positive Integer scalar)
%           
%           ('MaxEdgeGap' -> Positive integer scalar)
% 
%           ('PreEdgeGrow' -> Positive integer scalar)
%
%       ('BatchMode' -> True/False)
%       If this option value is set to true, all graphical output and user
%       interaction is suppressed.
%
% 
% Output:
% 
%   movieData - The modified movieData structure, with the mask refinement
%   logged in it, including all parameters used.
% 
%   The refined masks will replace the old masks in each channel's mask
%   directory.
% 
% 
% Hunter Elliott 
% 1/2010
%

%% ----------- Input --------- %%

movieData = setupMovieData(movieData,'maskRefinement');

[batchMode,iChannels,doCleanUp,doEdgeRefine,...
    minSize,closeRad,nObjects,...
    maxEdgeAdjust,maxEdgeGap,preEdgeGrow,fillHoles] = parseInput(varargin);


% ----- Default Values ----- %

if isempty(batchMode)
    batchMode = false;
end

if isempty(iChannels)
    iChannels = selectMovieChannels(movieData,1,'Select the channel(s) for mask refinement:');
end

if isempty(doCleanUp)
    doCleanUp = true;
end

if isempty(doEdgeRefine)
    doEdgeRefine = false;
end

if isempty(minSize)
    minSize = 10;
end

if isempty(closeRad)
    closeRad = [0 3];
end

if isempty(nObjects)
    nObjects = 1;
end

if isempty(fillHoles)
    fillHoles = true;
end

if ~checkMovieMasks(movieData,iChannels)
    error('Selected channels must have valid masks! Check channel masks and selected channels!')
end

if ~(doCleanUp || doEdgeRefine)
    error('You must enable mask cleanup and/or edge refinement! Otherwise this function has nothing to do!')
end


%% ----------- Init ----------- %%


maskNames = getMovieMaskFileNames(movieData,iChannels);

if doEdgeRefine %Images are only needed for edge-refinement
    imageNames = getMovieImageFileNames(movieData,iChannels);
end

nChan = length(iChannels);

if closeRad(1) > 0 %If closure is to be performed, create the structuring element
    seClose1 = strel('disk',closeRad(1));    
end

if closeRad(2) > 0
    seClose2 = strel('disk',closeRad(2));
end


%% ---------------- Mask Refinement --------------- %%


disp('Starting mask refinement...')

for iChan = 1:nChan

    nImages = movieData.nImages(iChannels(iChan));
    
    disp(['Refining masks for channel ' movieData.channelDirectory{iChannels(iChan)}]);
    if ~batchMode
        wtBar = waitbar(0,['Please wait, refining masks for channel '...
            movieData.channelDirectory{iChannels(iChan)} '....']);
    end

    
    for iImage = 1:nImages;
        
        %Load the mask for this frame/channel
        currMask = imread(maskNames{iChan}{iImage});
        
        
        % ----- Mask Clean Up ------ %
        
        if doCleanUp
            
            %Perform initial closure operation
            if closeRad(1) > 0
                currMask = imclose(currMask,seClose1);            
            end
            
            %Remove objects that are too small
            if ~isinf(minSize)                
                currMask = bwareaopen(currMask,minSize);                                       
            end
        
            %Connect fractured objects
            if closeRad(2) > 0
                currMask = imclose(currMask,seClose2);                                
            end
            
           
        end
        
        
        % --------- Mask Edge-Refinement ------ %
        if doEdgeRefine
            
            %Load the current image
            currImage = imread(imageNames{iChan}{iImage});
            
            %Call the edge-refinement function
            currMask = refineMaskEdges(currMask,currImage,...
                maxEdgeAdjust,maxEdgeGap,preEdgeGrow);
            
            
        end
        % ---------- Object Selection -------- %
        
        %Keep only the largest objects
        if doCleanUp && ~isinf(nObjects)
                
            %Label all objects in the mask
            labelMask = bwlabel(currMask);

            %Get their area
            obAreas = regionprops(labelMask,'Area');      

            %First, check that there are objects to remove
            if length(obAreas) > nObjects 
                obAreas = [obAreas.Area];
                %Sort by area
                [~,iSort] = sort(obAreas,'descend');
                %Keep only the largest requested number
                currMask = false(size(currMask));
                for i = 1:nObjects
                    currMask = currMask | labelMask == iSort(i);
                end
            end
        end
        
        % ------ Hole-Filling ----- %
        if fillHoles
            currMask = imfill(currMask,'holes');
        end
        
        %Write the refined mask to file, over-writing the previous mask.
        imwrite(currMask,maskNames{iChan}{iImage});
        
        if ~batchMode
            waitbar(iImage/nImages,wtBar)
        end        
        
    end
    
    if ~batchMode && ishandle(wtBar)
        close(wtBar);
    end    
end

%% ------ Finalize ------ %%


%Store parameters/settings in movieData structure
movieData.maskRefinement.status = 1;
movieData.maskRefinement.dateTime = datestr(now);
movieData.maskRefinement.doCleanUp = doCleanUp;
movieData.maskRefinement.doEdgeRefine = doEdgeRefine;
movieData.maskRefinement.minSize = minSize;
movieData.maskRefinement.closeRad = closeRad;
movieData.maskRefinement.nObjects = nObjects;
movieData.maskRefinement.maxEdgeAdjust = maxEdgeAdjust;
movieData.maskRefinement.maxEdgeGap = maxEdgeGap;
movieData.maskRefinement.preEdgeGrow = preEdgeGrow;
movieData.maskRefinement.iFrom(iChannels) = iChannels;

%Save the updated movieData structure.
updateMovieData(movieData);

if ~batchMode && ishandle(wtBar) %In case the user closed it
    close(wtBar)
end


function [batchMode,iChannels,doCleanUp,doEdgeRefine,minSize,closeRad, ...
         nObjects,maxEdgeAdjust,maxEdgeGap,preEdgeGrow,fillHoles] = parseInput(argArray)


%Init output

batchMode = [];
iChannels = [];
doCleanUp = [];
doEdgeRefine = [];
minSize = [];
closeRad = [];
nObjects = [];
maxEdgeAdjust =[];
maxEdgeGap = [];
preEdgeGrow = [];
fillHoles = [];

if isempty(argArray)
    return
end

nArg = length(argArray);

%Make sure there is an even number of arguments corresponding to
%optionName/value pairs
if mod(nArg,2) ~= 0
    error('Inputs must be as optionName / value pairs!')
end

for i = 1:2:nArg
    
   switch argArray{i}                     
              
       case 'BatchMode'
           batchMode = argArray{i+1};
           
       case 'ChannelIndex'
           iChannels = argArray{i+1};
           
       case 'MaskCleanUp'
           doCleanUp = argArray{i+1};
           
       case 'EdgeRefinement'
           doEdgeRefine = argArray{i+1};
           
       case 'MinimumSize'
           minSize = argArray{i+1};
           
       case 'ClosureRadius'
           closeRad = argArray{i+1};
           
       case 'ObjectNumber'
           nObjects = argArray{i+1};
           
       case 'MaxEdgeAdjust'
           maxEdgeAdjust = argArray{i+1};
           
       case 'MaxEdgeGap'
           maxEdgeGap = argArray{i+1};
           
       case 'PreEdgeGrow'
           
           preEdgeGrow = argArray{i+1};
       
       case 'FillHoles'
           
           fillHoles = argArray{i+1};
           
           
       otherwise
       
           error(['"' argArray{i} '" is not a valid option name! Please check input!'])
   end
   
end