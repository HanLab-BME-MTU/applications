classdef CellInlierDetectionClassifier_OneClassSVM < handle
    
   properties(SetAccess = private)
      wekaModel; 
   end
    
   methods
       
       function obj = CellInlierDetectionClassifier_OneClassSVM( wekaModelPath )
           
           AddWekaClassesToPath( {'libsvm'} );
           
           % weka imports
           import weka.core.*;    

           obj.wekaModel = weka.core.SerializationHelper.readAll( wekaModelPath );
           
       end
       
       function [predictedLabel, predictionProbability] = predictCell(obj, imInput, imLabelCellSeg, cellId, spacing)
           
            curFeatureStruct = ComputeCellInlierDetectionFeatures( imInput, imLabelCellSeg, cellId, spacing );
            [featureVec , featureNameList] = ConvertFeatureStructToFeatureVec( curFeatureStruct );        

            [predictedLabel, predictionProbability] = obj.predictInstance(featureVec, featureNameList);
           
       end
           
   end
   
   methods (Access = private)

       function [predictedLabel, predictionProbability] = predictInstance(obj, featureVec, featureNameList)
           
          [wekaInstance, wekaAttributeList, wekaModelClassNameList] = ConvertFeatureVectorToWekaInstance(obj.wekaModel, featureVec, featureNameList);
           
          curDataSet = weka.core.Instances('TestInstance', wekaAttributeList, 1 );
          curDataSet.add(wekaInstance);
          curDataSet.setClassIndex(0);
          
          wekaModelClassifier = obj.wekaModel(1);
          wekaModelAttributeHeader = obj.wekaModel(2);
          
          if wekaModelAttributeHeader.numAttributes ~= curDataSet.numAttributes || ~wekaModelAttributeHeader.equalHeaders(curDataSet)
              error( 'ERROR: mismatch between computed features and classifier features' );
          end
          
          predictedLabel = wekaModelClassifier.distributionForInstance(curDataSet.instance(0));
          predictionProbability = 1.0;
          
      end
       
end