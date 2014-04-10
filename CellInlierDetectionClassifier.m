classdef CellInlierDetectionClassifier < handle
    
   methods(Abstract)
       
       [predictedLabel, predictionProbability] = predictCell(obj, imageData, imLabelCellSeg, cellId, spacing);
       
   end
   
   methods(Static)
       
       function [imageDataAdjusted] = preprocessImageData(imageData)

           % standardize
           imageDataAdjusted = mat2gray(imageData) * 4096;

           % denoise
           imageDataAdjusted = matitk( 'FMEDIAN', [1,1,1], imageDataAdjusted );    

       end       

   end

end