function [ output_args ] = quickFixForControl052( analInfo,protS )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

frames2Fix  = find(arrayfun(@(x) ~isempty(analInfo(x).outlierFixFlag),1:length(analInfo))); 
for i = 1:length(frames2Fix)
   frameC =  frames2Fix(i) ; 
   fullMask = analInfo(frameC).masks.neuriteEdge;
   roiYXAll = bwboundaries(fullMask);
    pixIndicesAll = sub2ind(size(fullMask),roiYXAll{1}(:,1),roiYXAll{1}(:,2));
    
    

    % load the old thin body pixels
    pixIndicesThinBody =  analInfo(frameC).bodyEst.pixIndThinBody;
    % anything for now that is not thin body is new thick body
    pixIndicesThickBodyNew = setdiff(pixIndicesAll,pixIndicesThinBody);
    
    thickBodyMask = zeros(size(fullMask));
    thickBodyMask(pixIndicesThickBodyNew)=1;
    thickBodyMask = logical(thickBodyMask);
    % record
    analInfo(frameC).masks.thickBodyMask = thickBodyMask;
    analInfo(frameC).bodyEst.pixIndThickBody = pixIndicesThickBodyNew;
    
end 
save([pwd filesep 'analInfoTestSave.mat'],'analInfo','-v7.3'); 

end

