function [ filoInfo ] = GCAextractFilopodiaFromSubRegions(filoInfo,subRoiMask)
%
% INPUT: 
% filoInfo: Nx1 structure output from the reconstruction for each frame 
%           and stored in analInfo, where N is the number of filopodia 
%           detected per frame. 
% subRoiMask: an NxM logical mask of the subRegion from which one would
%            like to extract the filopodia 
% 
% 
% OUTPUT: a truncated filoInfo structure containing only those filopodia 
%        spatially restricted within that subRoiMask
% 
% CHECK INPUT: 



%% START
% first filter the filopodia for those attached to the body. 
% the field .type indicates this (0 is attached directly to the body)
 
filoType = vertcat(filoInfo(:).type); % extract types 
filoInfo = filoInfo(filoType==0); % filter by body

% get all attachment point coordinates of all filopodia in pixel idx
% (remember that you stored it such that the base is the first coordinate. 
  idxAttachFilo = arrayfun(@(x) filoInfo(x).Ext_pixIndices(1),1:length(filoInfo)); 

 % the way I set this up currently might have been a bit stupid.. as I 
 % don't know if the attachment points are actually on the mask or 1 pixel 
 % removed.  
 % check here. 
 subRoiMask = imdilate(subRoiMask,strel('disk',1)); 
 IN = subRoiMask(idxAttachFilo); % logical 1 if in region 0 is not. 
 
 filoInfo = filoInfo(IN); % get the filoInfo for the subregion 
 p.plots = 1; 
 if p.plots == 1 
 imshow(roiMask,[]); 
 hold on 
 arrayfun(@(x) plot(x.xyCoords(:,1),x.xyCoords(:,2),'color','r'),filoInfo); 
 end 
 
 
 
 
 
 
  

end

