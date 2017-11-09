function [erodMask,saveMask] = gcaMorphologicalOpeningWithGeometryConstraints(mask,varargin)
% 
% Small function that will mark which CCs from the erosion can be removed
% and NOT break the cell body these will be marked as potential higher
% confidence linkages.
% the CCs that should not be removed because they will create major
% breakages in the overal body ( could likely here filter by size as well -
% if a very small linkage or very low intensity linkage can remove)
% was findBreakageCCsFromErosion.m until renamed 20150511 and fixed input
%%  INPUT: 
% mask (REQUIRED) RxC logical array
%       of the original local thresholded image before morphological opening
%       : where R is the height (ny) and C is the width (nx) of the input image
%
% %% PARAMS: Morphological Opening %%
%
%    'DiskSizeLarge' (PARAM) : Positive scalar 
%       
%       Default = 6 for LifeAct: 4 for membrane markers
%       This parameter specifies the radius of the disk (in pixels) to be used for
%       the removal of thin objects from the input mask (ie the filopodia)
%       Larger values will remove thicker structures. Note for the lifeAct
%       channels that have very strong filopodia signal very often these
%       gradients for crossing filopodia tend not to be well segmented so
%       disks use are typically slightly bigger than if using a membrane
%       marker where the filpodia often exhibit weak signal relative to
%       entire image and are those not segmented at all.
%
%
%    'DiskSizeSmall' (PARAM) : Positive scalar 
%     
%     Default = 3 : 
%     Morphological 
%        
% 
%
%
%
%OUTPUT:
% erodMask: 
% 
% saveMask: 
%% INPUTPARSER
%%Input check
ip = inputParser;

ip.CaseSensitive = false;
ip.KeepUnmatched = true; 
%REQUIRED
ip.addRequired('mask');

% PARAMETERS
%ip.addParameter('TSOverlays',true,@(x) islogical(x));


ip.addParameter('DiskSizeLarge',6,@(x) isscalar(x));
ip.addParameter('DiskSizeSmall',3,@(x) isscalar(x)); 

ip.parse(mask,varargin{:});
 
 
%%
% my personal sanity check : decide if want to make part of trouble shoot
sanityCheck = 0;
%% Perform morphological opening
erodForBody = imopen(mask,strel('disk',ip.Results.DiskSizeLarge,0));

diffMask = mask-erodForBody;

CCDiff = bwconncomp(diffMask);
CCStart = bwconncomp(mask);
numStart = CCStart.NumObjects;

%% Check local morphological opening will result in an increase in the number of total connected components
% if so reduce the morphological opening radius for these regions. 
for iPiece = 1:CCDiff.NumObjects
    [ny,nx] = size(mask);
    maskPiece = zeros(ny,nx);
    
    maskPiece(CCDiff.PixelIdxList{iPiece}) = 1;
    
    
    
    CCEnd = bwconncomp(mask-maskPiece);
    numEnd= CCEnd.NumObjects;
    delta =  numStart - numEnd;
   
    if delta < 0 ; % the number of end pieces  more than the number of start pieces
        CCDiff.diffMark{iPiece}= 1;
    else
        CCDiff.diffMark{iPiece}= 0;
    end
%% Personal Sanity Check Plot: Not sure I want to make it an option yet
% as it is for EACH Piece so it makes a ton of figures
    % if make plot
    if sanityCheck == 1
        fsFigure(0.75)
        subplot(1,2,1);
        % make a label diff Mask
        spy(mask,'y')
        title(num2str(numStart));
        hold on
        spy(maskPiece,'r');
        subplot(1,2,2);
        spy(mask-maskPiece,'b');
        title(num2str(numEnd));
        saveas(gcf,[num2str((iPiece),'%03d') '.tif']);
        saveas(gcf,[num2str((iPiece),'%03d') '.fig']);
        close gcf
    end
    
end % for iPiece

%% Add back the connected component regions marked above  
filter = vertcat(CCDiff.diffMark{:});
filter = logical(filter);
saveMask = zeros(ny,nx);

% erosion pieces to save
saveMask(vertcat(CCDiff.PixelIdxList{filter})) = 1;

erodMask = (erodForBody|saveMask);
%% perfrom a smaller radius morphological opening to remove any filopodia like-structures along these pieces
erodMask = imopen(erodMask,strel('disk',ip.Results.DiskSizeSmall,0));
end

