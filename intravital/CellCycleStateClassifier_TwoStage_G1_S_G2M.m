classdef CellCycleStateClassifier_TwoStage_G1_S_G2M < CellCycleStateClassifier
    
   properties(SetAccess = private)
        wekaModel_Stage1;
        wekaModel_Stage2;
   end
    
   methods
       
       function obj = CellCycleStateClassifier_TwoStage_G1_S_G2M( wekaModelPath_Stage1, wekaModelPath_Stage2 )
           
           AddWekaClassesToPath();
           
           % weka imports
           import weka.core.*;    

           obj.wekaModel_Stage1 = weka.core.SerializationHelper.readAll( wekaModelPath_Stage1 );
           obj.wekaModel_Stage2 = weka.core.SerializationHelper.readAll( wekaModelPath_Stage2 );
           
       end
       
       function [predictedLabel, predictionProbability] = predictCell(obj, imageData, imValidROIMask, imLabelCellSeg, cellId, spacing)
           
            curFeatureStruct = ComputeCellStateClassificationFeatures( imageData, imValidROIMask, imLabelCellSeg, cellId, spacing );
            [featureVec , featureNameList] = ConvertFeatureStructToFeatureVec( curFeatureStruct );        

            % stage-1 
            [predictedLabel, predictionProbability] = obj.predictInstance(obj.wekaModel_Stage1, featureVec, featureNameList);
           
            % stage-2
            if ~ismember( predictedLabel, { 'G1', 'S', 'G2-M' } )
                error( 'Invalid class label - %s - received out of stage-1 classifier. Must be either G1 or S or G2-M', predictedLabel );
            end
            
            if strcmp( predictedLabel, 'G2-M' )
                [predictedLabel, predictionProbability] = obj.predictInstance(obj.wekaModel_Stage2, featureVec, featureNameList);
            end
       end
       
   end
   
   methods (Access = private)

       function [predictedLabel, predictionProbability] = predictInstance(obj, wekaModel, featureVec, featureNameList)
           
          [wekaInstance, wekaAttributeList, wekaModelClassNameList] = ConvertFeatureVectorToWekaInstance(wekaModel, featureVec, featureNameList);
           
          curDataSet = weka.core.Instances('TestInstance', wekaAttributeList, 1 );
          curDataSet.add(wekaInstance);
          curDataSet.setClassIndex(0);
          
          wekaModelClassifier = wekaModel(1);
          wekaModelAttributeHeader = wekaModel(2);
          
          if ~wekaModelAttributeHeader.equalHeaders(curDataSet)
              error( 'ERROR: mismatch between computed features and classifier features' );
          end
          
          wekaPredictedClassLabel = wekaModelClassifier.classifyInstance(curDataSet.instance(0)) + 1;
          wekaPredictionProbabilityDist = wekaModelClassifier.distributionForInstance( curDataSet.instance(0) );
          
          predictedLabel = wekaModelClassNameList{wekaPredictedClassLabel};
          predictionProbability = wekaPredictionProbabilityDist(wekaPredictedClassLabel);
          
       end
       
   end
   
end