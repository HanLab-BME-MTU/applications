classdef CellCycleStateClassifier_TwoStage_G1_SG2M < CellCycleStateClassifier
    
   properties(SetAccess = private)
        wekaModel_Stage1;
        wekaModel_Stage2;
        classNameList;
        classNameList_Stage1;
        classNameList_Stage2;
   end
    
   methods
       
       function obj = CellCycleStateClassifier_TwoStage_G1_SG2M( wekaModelPath_Stage1, wekaModelPath_Stage2 )
           
           AddWekaClassesToPath();
           
           % weka imports
           import weka.core.*;    

           obj.wekaModel_Stage1 = weka.core.SerializationHelper.readAll( wekaModelPath_Stage1 );
           obj.wekaModel_Stage2 = weka.core.SerializationHelper.readAll( wekaModelPath_Stage2 );
          
           classNameList = { 'G1', 'S', 'G2', 'M' };           
           
       end
       
       function [predictedLabel, classPredictionProbabilities] = predictCell(obj, imageData, imValidROIMask, imLabelCellSeg, cellId, spacing)
           
            curFeatureStruct = ComputeCellStateClassificationFeatures( imageData, imValidROIMask, imLabelCellSeg, cellId, spacing );
            [featureVec , featureNameList] = ConvertFeatureStructToFeatureVec( curFeatureStruct );        

            % stage-1 
            [predictedLabel, classProbs_Stage1] = obj.predictInstance(obj.wekaModel_Stage1, featureVec, featureNameList);
           
            % stage-2
            if ~ismember( predictedLabel, { 'G1', 'S-G2-M' } )
                error( 'Invalid class label - %s - received out of stage-1 classifier. Must be either G1 or S-G2-M', predictedLabel );
            end

            if strcmp(predictedLabel, 'S-G2-M')
                [predictedLabel, classProbs_Stage2] = obj.predictInstance(obj.wekaModel_Stage2, featureVec, featureNameList);
            end
            
       end
       
       function classNameList = getClassNameList(obj)
            classNameList = obj.classNameList;
       end
       
   end
   
   methods (Access = private)

       function [predictedLabel, classPredictionProbabilities] = predictInstance(obj, wekaModel, featureVec, featureNameList)
           
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
          predictedLabel = wekaModelClassNameList{wekaPredictedClassLabel};
          
          classPredictionProbabilities = wekaModelClassifier.distributionForInstance( curDataSet.instance(0) );
          
       end
       
   end
   
end