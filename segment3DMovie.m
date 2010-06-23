function movieData = segment3DMovie(movieData,paramsIn)
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
%   movieData - MovieData3D object describing the movie to segment.
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
% 
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
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
%               'HuntThresh' - Global threshold based on my histogram
%               analysis function thresholdFluorescenceImage.m
% 
%               'Gradient' - The image gradient is calculated
%               and this gradient image is then thresholded using
%               thresholdFluorescenceImage.m
%
% 
%       ('p.BatchMode' -> logical) If true, all graphical output is suppressed. 
%       Default is false.       
%
%       ('p.PostProcess' -> logical) If true, post-processing will be
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
%           ('NumObjects' -> positive integer) Keep only the largest p.NumObjects
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
gSig = 1; %Sigma of the filter used in the smoothed gradient filter.
dName = 'masks_channel_'; %Name for mask directories

%% ------ Input ----- %%


if ~isa(movieData,'MovieData3D')
    error('The first input argument must be a valid MovieData object!')
end

%Check for existing seg processes
iSegProc = cellfun(@(x)(isa(x,'SegmentationProcess3D')),movieData.processes_);

if isempty(iSegProc)
    iSegProc = numel(movieData.processes_)+1;
    movieData.addProcess(SegmentationProcess3D(movieData));
end

if nargin < 2
    paramsIn = [];
end

p = parseProcessParams(movieData.processes_{iSegProc},paramsIn);

if isempty(p.ChannelIndex)
    %Default is to use all channels
    p.ChannelIndex = 1:length(movieData.channelDirectory);
end


%% ----- Init -----%%

%Get number of channels to segment
nChanSeg = length(p.ChannelIndex);

%Set up mask output paths
for j = 1:nChanSeg
    movieData.processes_{iSegProc}.setMaskPath(p.ChannelIndex(j),...
            [p.OutputDirectory filesep dName num2str(p.ChannelIndex(j))]);
        
    mkClrDir(movieData.processes_{iSegProc}.maskPaths_{p.ChannelIndex(j)});
end

%Create structuring element for post-processing
if p.PostProcess && p.ClosureRadius > 0
    closeBall = binarySphere(p.ClosureRadius);
end



%% ---------- Segmentation ----------%%
%Go through each image in each channel and create a mask

if ~p.BatchMode
    wtBar = waitbar(0,'Please wait, creating masks....');
end

nImages = movieData.nFrames_;
nImTot = nImages * nChanSeg;

imNames = movieData.getImageFileNames(p.ChannelIndex);

for iChan = 1:nChanSeg                
    
    maskDir = movieData.processes_{iSegProc}.maskPaths_{p.ChannelIndex(iChan)};
    imDir = movieData.channelPath_{p.ChannelIndex(iChan)};
    
    if p.FixJumps
        threshVals = nan(nImages,1);                                
    end
    
    
    for iImage = 1:nImages

        %Load the current image
        currIm = stackRead([imDir filesep imNames{iChan}{iImage}]);
               
        % ---- Perform initial segmentation ---- %
        
        switch p.Method
                    
            case 'Otsu'
                
                %Perform otsu thresholding to get the mask
                currThresh = graythresh(currIm);
                range = getrangefromclass(currIm);
                currThresh = range(2) * currThresh; %Convert the stupid fractional threshold
                
            case 'HuntThresh'
                
                try                
                    currThresh = thresholdFluorescenceImage(currIm);
                    
                 catch errMess %If the auto-thresholding fails, 
                    if p.FixJumps
                        % just force use of last good threshold, if Fix
                        % Jumps enabled
                        currThresh = Inf;                        
                    else
                        %Otherwise, throw an error.
                        error(['Error: ' errMess.message ' Try enabling the FixJumps option!'])                                                
                    end                    
                end
                
            case 'Gradient'
                
                %Get the gradient of the image
%                 [gX,gY,gZ] = gradient(double(currIm));
%                 currIm = sqrt( gX .^2 +  gY .^2 + gZ .^2); %Just overwrite the image, save memory etc.                              
%                 
                currIm = matitk('FGMS',gSig,double(currIm));
                
                %Threshold this gradient based on intensity histogram
                try
                    currThresh = thresholdFluorescenceImage(currIm);                
                catch errMess %If the auto-thresholding fails, 
                    if p.FixJumps
                        % just force use of last good threshold, if Fix
                        % Jumps enabled
                        currThresh = Inf;                        
                    else
                        %Otherwise, throw an error.
                        error(['Error: ' errMess.message ' Try enabling the FixJumps option!'])                                                
                    end                    
                end
                
            otherwise
                
                error(['The segmentation method ' p.Method ' is not recognized! Please check Method option!'])
                                                
        end     
                
        
        
        % --- Check the new thresholdValue if requested ---- % 
        if p.FixJumps
            
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
        
        if p.PostProcess                                
            
            currMask = bwareaopen(currMask,p.MinVolume);
            
            currMask = imclose(currMask,closeBall);
                    
            labelMask = bwlabeln(currMask);
            
            if p.NumObjects > 0
                rProp = regionprops(labelMask,'Area');
                [goodObj,iGood] = sort([rProp(:).Area],'descend'); %#ok<ASGLU>

                currMask = false(size(currMask));
                for i = 1:p.NumObjects
                    currMask = currMask | (labelMask == iGood(i));
                end
            end
        end
           
        %We want to compress the masks, so don't use stackWrite.m
        for i = 1:size(currIm,3)
            %Append each z-slice to the tiff
            imwrite(currMask(:,:,i),[maskDir filesep ...
            'mask_' imNames{iChan}{iImage}(1:end-3) 'tif'],'tif','WriteMode','append')
        end
        
        if ~p.BatchMode
            waitbar( (iImage+(nImages*(iChan-1))) / nImTot,wtBar)
        end
    end
end

%% ------Output/Finalization----- %%

movieData.processes_{iSegProc}.setSuccess(true);
movieData.processes_{iSegProc}.setDateTime;
movieData.saveMovieData; %Save the new movieData to disk

disp('Finished thresholding!')

if ~p.BatchMode && ishandle(wtBar)
    close(wtBar)
end

