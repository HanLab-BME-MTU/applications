function maskProp = adjustROIMaskGeometry(maskProp,roiInf,pixXY,pixZ)
%Adjusts the xyz coordinates to agree with the new origin of the ROI.

newOrig = [roiInf.cropX(1) roiInf.cropY(1) roiInf.cropZ(1)]-1;%Offset used to adjust for new image origin location
newOrig(3) = newOrig(3) * pixZ/pixXY;   %Adjust z for symmetric voxels

adjFun = @(x)(x - repmat(newOrig([2 1 3]),[size(x,1),1]));

maskProp.SmoothedSurface.vertices = adjFun(maskProp.SmoothedSurface.vertices);
maskProp.CenterMost = adjFun(maskProp.CenterMost);
maskProp.Centroid = adjFun(maskProp.Centroid);


