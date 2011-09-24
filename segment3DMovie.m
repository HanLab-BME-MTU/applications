function movieData = segment3DMovie(movieData,paramsIn)
%SEGMENT3DMOVIE creates 3D masks for the input 3D movie
% 
% movieData = segment3DMovie(movieData)
% 
% movieData = segment3DMovie(movieData,paramsIn)
% 
% This function creates 3D masks for the input movie. The masks can be
% initially created using a few methods (see below). Post-processing is
% also performed if requested. This function is designed to segment
% fluorescence images, using either a volumetric or membrane-bound label.
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
%       ('ThresholdValue'->Positive scalar) If specified, all masks will be
%       created by thresholding at this value. If not specified, one of the
%       automatic methods will be used. See 'Method' option below.
%
%       ('Method' - Character array) The name of the method to use for
%       segmentation. 
%
%           Available methods are:
%
%               'Otsu' - Global Otsu thresholding.
%
%               'HuntThresh' - Global threshold based on my histogram
%               analysis function thresholdFluorescenceImage.m This is the
%               default.
% 
%               'Gradient' - The image gradient is calculated
%               and this gradient image is then thresholded using
%               thresholdFluorescenceImage.m
%
%       ('PreFilterSig'->positive scalar) Specifies the sigma of the
%       gaussian filter to apply the images prior to segmentation.
%       Optional. Default is 0 (no filtering).
% 
%       ('BatchMode' -> logical) If true, all graphical output is suppressed. 
%       Default is false.       
%
%       ('PostProcess' -> logical) If true, post-processing will be
%       performed on the mask. This includes closure, area-opening and
%       object selection. Default is true.
%
%       ('FixJumps' -> Logical) If true, sudden changes in the threshold
%       values used will be suppressed by using the last good threshold
%       value.
%       Default is true. 
%
%
%       Post Processing Option/Value pairs:
%            
%           ('MinVolume' -> Positive integer) The minimum volume (in
%           voxels) of objects to keep. Objects smaller than this will be
%           removed from masks. If zero, no objects are removed. Default is
%           25 voxels.
% 
%           ('NumObjects' -> positive integer) Keep only the largest p.NumObjects
%           objects in the masks. If set to zero, all objects are kept.
%           Default is 1.
%
%           ('ClosureRadius' -> positive integer) Radius of structuring
%           element to use in closure operation. If zero, no closure is
%           performed. Default is 2 pixels
%
%           ('FillHoles' -> 0,1,2,3 or 4) Specifies if and how to do hole-filling
%           in the mask:
%               0 - No hole filling.
%               1-3 - Lenient hole filling. Does hole-filling in 2D along
%               each dimension of the mask, and then any pixels which were
%               filled in the specified number of dimensions. That is, if
%               FillHoles = 2, then only pixels which were filled in 2 or
%               more of the 2D hole-fills will be filled in the final mask.
%               A closure operation will then be applied to this filled mask.
%               4 - Strict hole filling. Only true 3D holes will be filled
%               (pixels which cannot be reached by filling in from the
%               image border)
%
% Hunter Elliott
% 11/2009
%

%% ----- Parameters ---- %%

maxJump = .25; %The maximum fractional change in a threshold value to allow if the FixJumps option is enabled. %TEMP - allow user-specification
%gSig = 1; %Sigma of the filter used in the smoothed gradient filter.
dName = 'masks_channel_'; %Name for mask directories


%% ------ Input ----- %%


if ~isa(movieData,'MovieData3D')
    error('The first input argument must be a valid MovieData object for a 3D movie!')
end

%Check for existing seg processes
iSegProc = movieData.getProcessIndex('SegmentationProcess3D',1,0);

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

%MORE INPUT CHECKING!?!


%% ----- Init -----%%

%Get number of channels to segment
nChanSeg = length(p.ChannelIndex);

%Set up mask output paths
maskDir = cell(1,nChanSeg);
for j = 1:nChanSeg
    
    maskDir{j} = [p.OutputDirectory filesep dName num2str(p.ChannelIndex(j))];
    movieData.processes_{iSegProc}.setOutMaskPath(p.ChannelIndex(j),maskDir{j});        
    mkClrDir(maskDir{j});
    
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
imDir = movieData.getChannelPaths(p.ChannelIndex);
    

for iChan = 1:nChanSeg                
    
    
    if p.FixJumps
        threshVals = nan(nImages,1);                                
    end
    
    
    for iImage = 1:nImages

        %Load the current image
        currIm = stackRead([imDir{iChan} filesep imNames{iChan}{iImage}]);
        
        
        %Pre-filter the image if requested
        if p.PreFilterSig > 0
            ogClass = class(currIm);
            currIm = filterGauss3D(double(currIm),p.PreFilterSig,'symmetric');
            currIm = cast(currIm,ogClass);%Cast back so Otsu will work...
        end
               
        % ---- Perform initial segmentation ---- %
        if isempty(p.ThresholdValue)
            switch p.Method

                case 'Otsu'

                    %Perform otsu thresholding to get the mask
                    currThresh = graythresh(currIm);
                    range = getrangefromclass(currIm);
                    currThresh = range(2) * currThresh; %Convert the stupid fractional threshold
                    
                    %Threshold the damn thing!
                    currMask = currIm > currThresh;

                case 'HuntThresh'

                    try                
                        currThresh = thresholdFluorescenceImage(currIm);

                     catch errMess %If the auto-thresholding fails, 
                        if p.FixJumps
                            % just force use of last good threshold, if Fix
                            % Jumps enabled
                            currThresh = Inf;                        
                        else
                            bsFill%Otherwise, throw an error.
                            error(['Error: ' errMess.message ' Try enabling the FixJumps option!'])                                                
                        end                    
                    end
                    
                    %Threshold the damn thing!
                    currMask = currIm > currThresh;        

                case 'Gradient'
                    %Get gradient of filtered image
                    [gX,gY,gZ] = gradient(double(currIm));
                    %and magnitude of gradient
                    currIm = sqrt( gX .^2 +  gY .^2 + gZ .^2);
                    
                    %There are often strong intensity variations in the
                    %first few stacks at top and bottom, so we just cut
                    %these out. - TEMP??! -HLE
                    nCut = 4;
                    currIm(:,:,1:nCut) = min(currIm(:));
                    currIm(:,:,end:-1:end-nCut) = min(currIm(:));                    
                        
                    %MATITK is slow, so I stopped using it - HLE
                    %Use MATITK for smoothed gradient calculation  
                    %currIm = matitk('FGMS',gSig,double(currIm));

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
                    
                    %Threshold the damn thing!
                    currMask = currIm > currThresh;        
                    
                case 'SurfaceEnhancement'
                    
                    %Surface-filtering and intensity segmentation.
                    currMask = huntersFancySegmentation3(currIm);

                otherwise

                    error(['The segmentation method ' p.Method ' is not recognized! Please check Method option!'])

            end                         
        

            % --- Check the new thresholdValue if requested ---- % 
            if p.FixJumps

                if iImage == 1
                    if currThresh == Inf
                        error('Failed on first frame: couldn'' automatically determine a threshold!');
                    else                    
                        threshVals(iImage) = currThresh; %Nothing to compare 1st frame to
                    end
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
                
        else
            currThresh = p.ThresholdValue;
        end
        
                
        % ---- Post-Processing if Requested ----- %
        
        if p.PostProcess     
            
            
            currMask = postProcess3DMask(currMask,p);
            
            
%             currMask = bwareaopen(currMask,p.MinVolume);
%             
%             currMask = imclose(currMask,closeBall);
%                     
%             labelMask = bwlabeln(currMask);
%             
%             if p.NumObjects > 0
%                 rProp = regionprops(labelMask,'Area');
%                 [goodObj,iGood] = sort([rProp(:).Area],'descend'); %#ok<ASGLU>
% 
%                 currMask = false(size(currMask));
%                 for i = 1:min(p.NumObjects,numel(iGood))
%                     currMask = currMask | (labelMask == iGood(i));
%                 end
%             end
%                                         
%             if p.FillHoles > 0 && p.FillHoles < 4
%                 
%                 %'Lenient' semi-2D hole filling                
%                 mFill1 = false(size(currMask));
%                 mFill2 = false(size(currMask));
%                 mFill3 = false(size(currMask));                
%                 for j = 1:size(currMask,1)
%                     mFill1(j,:,:) = imfill(squeeze(currMask(j,:,:)),'holes');                    
%                 end
%                 for j = 1:size(currMask,2)
%                     mFill2(:,j,:) = imfill(squeeze(currMask(:,j,:)),'holes');                    
%                 end                    
%                 for j = 1:size(currMask,3)
%                     mFill3(:,:,j) = imfill(currMask(:,:,j),'holes');                    
%                 end   
%                 %Only fill pixels which were filled in the specified number
%                 %of dimensions
%                 currMask = double(mFill1)+double(mFill2)+double(mFill3) >= p.FillHoles;                
%                 %Follow this with another round of closure
%                 currMask = imclose(currMask,closeBall);                
%                     
%             elseif p.FillHoles == 4
%                 %Do strict, full-3D hole-filling.
%                 currMask = imfill(currMask,'holes');                
%             end
        end
           
        %We want to compress the masks, so don't use stackWrite.m
        for i = 1:size(currIm,3)
            %Append each z-slice to the tiff. This uses matlab default
            %compression for binary files - bitpacking.
            imwrite(currMask(:,:,i),[maskDir{iChan} filesep ...
            'mask_' imNames{iChan}{iImage}(1:end-3) 'tif'],'tif','WriteMode','append')
        end
        
        if ~p.BatchMode
            waitbar( (iImage+(nImages*(iChan-1))) / nImTot,wtBar)
        end
    end
end

%% ------Output/Finalization----- %%


movieData.processes_{iSegProc}.setDateTime;
movieData.save; %Save the new movieData to disk

disp('Finished thresholding!')

if ~p.BatchMode && ishandle(wtBar)
    close(wtBar)
end

