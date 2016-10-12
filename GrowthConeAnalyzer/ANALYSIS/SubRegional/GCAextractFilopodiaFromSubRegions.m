function [ filoInfoSub] = GCAextractFilopodiaFromSubRegions(filoInfo,subRoiMask)
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
filoInfoBody = filoInfo(filoType==0 | filoType ==1); % filter by body

% get all attachment point coordinates of all filopodia in pixel idx
% (remember that you stored it such that the base is the first coordinate.
idxAttachFilo = arrayfun(@(x) filoInfoBody(x).Ext_pixIndices(1),1:length(filoInfoBody));

% the way I set this up currently might have been a bit stupid.. as I
% don't know if the attachment points are actually on the mask or 1 pixel
% removed.
% check here.
subRoiMask = imdilate(subRoiMask,strel('disk',2));
idxAttachFilo = idxAttachFilo(~isnan(idxAttachFilo));
IN = subRoiMask(idxAttachFilo); % logical 1 if in region 0 is not.

filoInfoSub = filoInfoBody(IN); % get the filoInfo for the subregion
typesInSub = vertcat(filoInfoSub(:).type);


if sum(typesInSub>0)
    % partition out the branch data
    filoInfoBranch = filoInfoSub(typesInSub>0);
    %filoInfoSub = filoInfoSub(typesInSub==0);
    
    
    % final all the filo in the subRoi that are of that group
    groupsToSave = vertcat(filoInfoBranch(:).groupCount);
    
    % for each branch group extract the data.
    for iBranch = 1:length(groupsToSave)
        idxBranchesFrames = arrayfun(@(x) filoInfo(x).groupCount==groupsToSave(iBranch) & filoInfo(x).type ~=1,1:length(filoInfo));
        idxBranches{iBranch} = find(horzcat(idxBranchesFrames));
        
    end
    
    idxBranches = horzcat(idxBranches{:});
    
    filoInfoSave = filoInfo(idxBranches);
    filoInfoSub = [ filoInfoSub  filoInfoSave]; 
end
% get the groupcount of any that are participating




p.plots = 0;
if p.plots == 1
    imshow(subRoiMask,[]);
    hold on
    %arrayfun(@(x) plot(x.Ext_coordsXY(:,1),x.Ext_coordsXY(:,2),'color','r'),filoInfo);
    arrayfun(@(x) plot(x.Ext_coordsXY(:,1),x.Ext_coordsXY(:,2),'color','g'),filoInfoSub);
end








end

