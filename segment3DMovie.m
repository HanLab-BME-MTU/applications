function movieData = segment3DMovie(movieData,varargin)
% 
% movieData = segment3DMovie(movieData)
% 
% movieData = segment3DMovie(movieData,'OptionName',optionValue)
% 
% This function creates 3D masks for the input movie. The masks can be
% initially created using a few methods (see below). Post-processing is also
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
%       (OptionName -> possibleValues)
%
%       ('ChannelIndex' -> positive integer) Integer indices of the
%       channel(s) to segment. Default is to segment all image channels.
%
%       ('Method' - Character array) The name of the method to use for
%       segmentation. 
%
%           Available methods are:
% 
%               'Otsu' - Global Otsu thresholding. This is the default.
% 
%               'Gradient' - The image gradient is calculated
%               and this gradient image is then thresholded using
%               thresholdFluorescenceImage.m
%
% 
%       ('BatchMode' -> logical) If true, all graphical output is suppressed. 
%       Default is false.       
%
%       ('PostProcess' -> logical) If true, post-processing will be
%       performed on the mask. This includes closure, area-opening and
%       object selection. Default is false.
%
%       ('FixJumps' -> Logical) If true, sudden changes in the threshold
%       values used will be suppressed by using the last good threshold
%       value.
%       Default is false. 
%
%
%       Post Processing Option/Value pairs:
%            
%           ('MinVolume' -> Positive integer) The minimum volume (in
%           voxels) of objects to keep. Objects smaller than this will be
%           removed from masks. If zero, no objects are removed. Default is
%           100 pixels.
% 
%           ('NumObjects' -> positive integer) Keep only the largest nObjects
%           objects in the masks. If set to zero, all objects are kept.
%           Default is 1.
%
%           ('ClosureRadius' -> positive integer) Radius of structuring
%           element to use in closure operation. If zero, no closure is
%           performed. Default is 3 pixels
%
% Hunter Elliott
% 11/2009
%

%% ----- Parameters ---- %%

maxJump = .25; %The maximum fractional change in a threshold value to allow if the FixJumps option is enabled


%% ------ Input ----- %%


movieData = setup3DMovieData(movieData,'masks');


[batchMode,iChannels,postProcess,minVolume,nObjects,closeRad,methName,fJump] = parseInput(varargin);



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

%Create structuring element for post-processing
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
    
    imNames = getMovieImageFileNames(movieData,iChannels);
    
    nImages = movieData.nImages(iChannels(iChan));
    
    %Check/create the directory for masks for this channel
    if ~isempty(movieData.masks.channelDirectory{iChan})
        maskDir = [movieData.masks.directory filesep movieData.masks.channelDirectory{iChan}];
    else
        maskDir = movieData.masks.directory; %This prevents addition of a trailing file seperator
    end
    
    if ~exist(maskDir,'dir')
        mkdir(maskDir)
    end
    
    nDig = floor(log10(nImages))+1;
    fString = ['%0' num2str(nDig) '.0f']; %Get formatting string for image names
    
    if fJump
        threshVals = nan(nImages,1);                                
    end
    
    
    for iImage = 1:nImages

        %Load the current image
        currIm = stackRead(imNames{iChan}{iImage});
               
        % ---- Perform initial segmentation ---- %
        
        switch methName
                    
            case 'Otsu'
                
                %Perform otsu thresholding to get the mask
                currThresh = graythresh(currIm);
                range = getrangefromclass(currIm);
                currThresh = range(2) * currThresh; %Convert the stupid fractional threshold                                           
                
            case 'Gradient'
                
                %Get the gradient of the image
                [gX,gY,gZ] = gradient(double(currIm));
                currIm = sqrt( gX .^2 +  gY .^2 + gZ .^2); %Just overwrite the image, save memory etc.                              
                
                %Threshold this gradient based on intensity histogram
                try
                    [currMask,currThresh] = thresholdFluorescenceImage(currIm);                
                catch errMess %If the auto-thresholding fails, 
                    if fJump
                        % just force use of last good threshold, if Fix
                        % Jumps enabled
                        currThresh = Inf;                        
                    else
                        %Otherwise, throw an error.
                        error(['Error: ' errMess.message ' Try enabling the FixJumps option!'])                                                
                    end                    
                end
                
            otherwise
                
                error(['The segmentation method ' methName ' is not recognized! Please check Method option!'])
                                                
        end     
                
        
        
        % --- Check the new thresholdValue if requested ---- % 
        if fJump
            
            if iImage == 1
                threshVals(iImage) = currThresh; %Nothing to compare 1st frame to
            else
                if abs(currThresh / threshVals(find(~isnan(threshVals),1,'last'))-1) > maxJump
                    %If the change was too large, don't store this threshold
                    %and use the most recent good value
                    threshVals(iImage) = NaN;
                    currThresh = threshVals(find(~isnan(threshVals),1,'last'));
                else
                    threshVals(iImage) = currThresh;
                end 
            end
        end
        
        
        %Threshold the damn thing!
        currMask = currIm > currThresh;
        
        
        
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
            'mask_' num2str(iImage,fString) '.tif'],'tif','WriteMode','overwrite')
        for i = 2:size(currIm,3)
            %Append the rest of the z-slices to this tiff
            imwrite(currMask(:,:,i),[maskDir filesep ...
            'mask_' num2str(iImage,fString) '.tif'],'tif','WriteMode','append')
        
        end
        
        if ~batchMode
            waitbar( (iImage+(nImages*(iChan-1))) / nImTot,wtBar) %it's approximate, but who cares
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

function [batchMode,iChannels,postProcess,minVolume,nObjects,closeRad,methName,fJump]...
    = parseInput(argArray)
%Sub-function for parsing input optionName/value pairs

%Default values
batchMode = false;
iChannels = [];
postProcess = false;
minVolume = 100;
nObjects = 1;
closeRad = 3;
methName = 'Otsu';
fJump = false;

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
              
       case 'BatchMode'
           batchMode = argArray{i+1};
           
       case 'ChannelIndex'
           iChannels = argArray{i+1};
           
       case 'PostProcess'
           postProcess = argArray{i+1};
           
       case 'MinVolume'
           minVolume = argArray{i+1};
           
       case 'NumObjects'
           nObjects = argArray{i+1};
           
       case 'CloseRad'
           nObjects = argArray{i+1};
           
       case 'Method'
           methName = argArray{i+1};
           
       case 'FixJumps'
           fJump = argArray{i+1};
           
       otherwise
       
           error(['"' argArray{i} '" is not a valid option name! Please check input!'])
   end
               
      
   
end





