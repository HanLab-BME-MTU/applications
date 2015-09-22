classdef CellCycleStateClassifier < handle
    
   methods(Abstract)
       
       [predictedLabel, classPredictionProbabilities, featureStruct] = predictCell(obj, imageData, imValidROIMask, imLabelCellSeg, cellId, spacing);
       [classNameList] = getClassNameList(obj);
       
   end
   
   methods(Static) 
       
       function imageDataAdjusted = preprocessImageData(imageData, flagParallelize)
          
           if ~exist( 'flagParallelize', 'var' )
                flagParallelize = false;
           end
           
           imageDataAdjusted = imageData;

           if flagParallelize
               
               parfor chid = 1:3

                   % standardize the intensity range
                   imageDataAdjusted{chid} = mat2gray( imageDataAdjusted{chid} ) * 4096;

                   % denoise using median filter
                   imageDataAdjusted{chid} = matitk( 'FMEDIAN', [1,1,1], imageDataAdjusted{chid} );    

               end
               
           else
               
               for chid = 1:3

                   % standardize the intensity range
                   imageDataAdjusted{chid} = mat2gray( imageDataAdjusted{chid} ) * 4096;

                   % denoise using median filter
                   imageDataAdjusted{chid} = matitk( 'FMEDIAN', [1,1,1], imageDataAdjusted{chid} );    

               end
               
           end
           
       end
           
   end
   
end