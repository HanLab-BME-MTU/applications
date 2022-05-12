function correctedField = correctLowFreqSigFromField(curField)
% function correctedField = correctLowFreqSigFromField(curField) background subtracts the low
% frequency signal in both x- and y- directions from the original field and
% overwrite the displacement field

% Look at x component of curField
[~,~,~,uxMap,uyMap] = generateHeatmapFromField(curField);
figure, imshow(uxMap,[])
figure, imshow(uyMap,[])


end

