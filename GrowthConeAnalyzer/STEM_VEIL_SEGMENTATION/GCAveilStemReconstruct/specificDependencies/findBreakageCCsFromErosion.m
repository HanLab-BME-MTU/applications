function [erodMask,saveMask] = findBreakageCCsFromErosion( mask,diskSize )
%Small function that will mark which CCs from the erosion can be removed
% and NOT break the cell body these will be marked as potential higher
% confidence linkages. 
% the CCs that should not be removed because they will create major
% breakages in the overal body ( could likely here filter by size as well -
% if a very small linkage or very low intensity linkage can remove)
%OUTPUT: 
sanityCheck = 0; 

 erodForBody = imopen(mask,strel('disk',diskSize,0));
 
 diffMask = mask-erodForBody; 
 
CCDiff = bwconncomp(diffMask);
CCStart = bwconncomp(mask); 
numStart = CCStart.NumObjects;
% this is likely computationally intensive but test 
for iPiece = 1:CCDiff.NumObjects
    [ny,nx] = size(mask); 
    maskPiece = zeros(ny,nx); 
    
    maskPiece(CCDiff.PixelIdxList{iPiece}) = 1; 
    
  
    
    CCEnd = bwconncomp(mask-maskPiece); 
    numEnd= CCEnd.NumObjects; 
   delta =  numStart - numEnd;
   
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
   
   
   
   
   if delta < 0 ; % the number of end pieces  more than the number of start pieces
       CCDiff.diffMark{iPiece}= 1; 
   else 
       CCDiff.diffMark{iPiece}= 0; 
   end 
end 
filter = vertcat(CCDiff.diffMark{:});
filter = logical(filter); 
saveMask = zeros(ny,nx); 

% erosion pieces to save
saveMask(vertcat(CCDiff.PixelIdxList{filter})) = 1; 
% Added 20140825 - helps to fill holes here (at least very small ones... so
% that don't reopen in the next step if there are just a few internal
% pixels) Example of the problem was with Frame 16 in Control 10 201212 
% ok doesn't work structures are still on the order of about 3 if get low
% to this size hard to keep connected - need to make rules or potentials
% with the intensity here was well. 
% filter out small erosion pieces 
%saveMask = imfill(saveMask,'holes'); % can be problematic really don't
%want to just blindly assume that the holes are a false negative often they
%can be due to overlapping filo that need to be eroded- better to lose them
%and then try to resave in the next step 
erodMask = (erodForBody|saveMask); 
 erodMask = imopen(erodMask,strel('disk',3,0)); 

end

