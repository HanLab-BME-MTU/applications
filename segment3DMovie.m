function movieData = segment3DMovie(movieData,varargin)
% 
% movieData = segment3DMovie(movieData)
% 
% movieData = segment3DMovie(movieData,'OptionName',optionValue)
% 
% This function creates 3D masks for the input movie. The masks are
% initially created using global Otsu thresholding. Post-processing is also
% performed if requested.
% 
% Input:
% 
%     movieData - Structure describing the movie, as created by
%     setup3DMovieData.m
% 
%     OptionName/value Pairs. - A string with option name followed by
%     its value.
%
%     Possible Option Names:
%
%       'channels' - Integer indices of the channel(s) to segment.
%       Default is to segment all image channels.
% 
% 
%       'batchMode' - If true, all graphical output is suppressed. 
%       Default is false.       
%
%       'postProcess' - If true, post-processing will be performed on the
%       mask. This includes closure, area-opening and object selection.
%       Default is false.
%
%       Post Processing Option/Value pairs:
%            
%           'minVolume' - The minimum volume (in voxels) of objects to
%           keep. Objects smaller than this will be removed from masks. If
%           zero, no objects are removed.
%           Default is 100 pixels.
% 
%           'nObjects' - Keep only the largest nObjects objects in the
%           masks. If set to zero, all objects are kept.
%           Default is 1.
%
%           'closureRadius' - Radius of structuring element to use in
%           closure operation. If zero, no closure is performed.
%           Default is 5 pixels
%
% Hunter Elliott, 11/2009
%

%% ------ Input ----- %%


movieData = setup3DMovieData(movieData,'masks');


[batchMode,iChannels,postProcess,minVolume,nObjects,closeRad] = parseInput(varargin);



if isempty(iChannels)
    %Default is to use all channels
    iChannels = 1:length(movieData.channelDirectory);
end


%% ----- Init -----%%

%Make sure this is a 3D movie
if ~isMovie3D(movieData)
    error('This function can only be used for 3D movies!')
end

%Get number of channels to segment
nChan = length(iChannels);

movieData.masks.channelDirectory(iChannels) = movieData.channelDirectory(iChannels);

movieData.masks.status = 0;

if postProcess && closeRad > 0
    closeDisk = strel('disk',closeRad);        
end

%% ---------- Segmentation ----------%%
%Go through each image in each channel and create a mask

if ~batchMode
    wtBar = waitbar(0,'Please wait, creating masks....');
end

nImTot = sum(movieData.nImages);


for iChan = 1:nChan

    imDir = [movieData.imageDirectory filesep movieData.channelDirectory{iChannels(iChan)}];
    imNames = dir([imDir filesep '*.STK']);

    nImages = length(imNames);
    
    %Check/create the directory for masks for this channel
    maskDir = [movieData.masks.directory filesep movieData.masks.channelDirectory{iChan}];
    if ~exist(maskDir,'dir')
        mkdir(maskDir)
    end

    for iImage = 1:nImages

        currIm = stkRead([imDir filesep imNames(iImage).name]);
        
        %Perform otsu thresholding to get the mask
        currMask = reshape(im2bw(currIm(:),graythresh(currIm)),size(currIm));
        
        % ---- Post-Processing if Requested ----- %
        if postProcess
        
            
            currMask = imclose(currMask,closeDisk);
            
            currMask = bwareaopen(currMask,minVolume);
                    
            labelMask = bwlabeln(currMask);
            
            rProp = regionprops(labelMask);
            [goodObj,iGood] = sort([rProp(:).Area],'descend');
            
            currMask = false(size(currMask));
            for i = 1:nObjects
                currMask = currMask | (labelMask == iGood(i));
            end
        end
        
        %Overwrite on the first z-slice to clear existing masks if present
        imwrite(currMask(:,:,1),[maskDir filesep ...
            'mask_' imNames(iImage).name(1:end-3) 'tif'],'WriteMode','overwrite')
        for i = 2:size(currIm,3)
            %Append the rest of the z-slices to this tiff
            imwrite(currMask(:,:,i),[maskDir filesep ...
            'mask_' imNames(iImage).name(1:end-3) 'tif'],'WriteMode','append')
        
        end
        
        if ~batchMode && mod(iImage,5) == 0
            waitbar( (iImage+(nImages*(iChan-1))) / nImTot,wtBar)
        end
    end
end

%% ------Output/Finalization----- %%

movieData.masks.status = 1;
movieData.masks.dateTime = datestr(now);


updateMovieData(movieData);

if ~batchMode && ishandle(wtBar)
    close(wtBar)
end

function [batchMode,iChannels,postProcess,minVolume,nObjects,closeRad]...
    = parseInput(argArray)
%Sub-function for parsing input optionName/value pairs

%Default values
batchMode = false;
iChannels = [];
postProcess = false;
minVolume = 100;
nObjects = 1;
closeRad = 5;


if isempty(argArray)
    return
end

nArg = length(argArray);

%Make sure there is an even number of arguments corresponding to
%optionName/value pairs
if mod(nArg,2) ~= 0
    error('Inputs must be as optionName/ value pairs!')
end

for i = 1:2:nArg
    
   switch argArray{i}                     
              
       case 'batchMode'
           batchMode = argArray{i+1};
           
       case 'channels'
           iChannels = argArray{i+1};
           
       case 'postProcess'
           postProcess = argArray{i+1};
           
       case 'minVolume'
           minVolume = argArray{i+1};
           
       case 'nObjects'
           nObjects = argArray{i+1};
           
       case 'closeRad'
           nObjects = argArray{i+1};
           
           
           
       otherwise
       
           error(['"' argArray{i} '" is not a valid option name! Please check input!'])
   end
               
      
   
end





