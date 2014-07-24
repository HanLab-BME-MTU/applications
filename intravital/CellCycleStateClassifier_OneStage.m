classdef CellCycleStateClassifier_OneStage < CellCycleStateClassifier
    
   properties(SetAccess = private)
        wekaModel;
   end
    
   methods
       
       function obj = CellCycleStateClassifier_OneStage( wekaModelPath )
           
           % weka imports
           import weka.core.*;    

           obj.wekaModel = weka.core.SerializationHelper.readAll( wekaModelPath );
           
       end
       
       function [predictedLabel, classPredictionProbabilities, featureStruct] = predictCell(obj, imageData, imValidROIMask, imLabelCellSeg, cellId, spacing)
           
            featureStruct = ComputeCellStateClassificationFeatures( imageData, imValidROIMask, imLabelCellSeg, cellId, spacing );
            [featureVec , featureNameList] = ConvertFeatureStructToFeatureVec( featureStruct );        

            [predictedLabel, classPredictionProbabilities] = wekaPredictInstance(obj.wekaModel, featureVec, featureNameList);
           
       end
       
       function classNameList = getClassNameList(obj)
           
           wekaModelAttributeHeader = obj.wekaModel(2);
           classNameList = cell(1, wekaModelAttributeHeader.attribute(0).numValues);
           for i = 1:numel(classNameList)
               classNameList{i} = char( wekaModelAttributeHeader.attribute(0).value(i-1) );
           end
           
       end
       
   end

end
