function [predictedLabel, classPredictionProbabilities] = wekaPredictInstance(wekaModel, featureVec, featureNameList)

    wekaModelClassifier = wekaModel(1);
    wekaModelAttributeHeader = wekaModel(2);
    
    [wekaInstance, wekaAttributeList, wekaModelClassNameList] = ConvertFeatureVectorToWekaInstance(wekaModel, featureVec, featureNameList);

    curDataSet = weka.core.Instances('TestInstance', wekaAttributeList, 1 );
    curDataSet.setClassIndex(0);

    [flagEqualHeaders, strErrorMessage] = wekaCheckIfHeadersAreEqual(wekaModelAttributeHeader, curDataSet);
    if ~flagEqualHeaders
      wekaModelAttributeHeader
      curDataSet
      error( 'ERROR: mismatch between computed features and classifier features\n%s\n', strErrorMessage );
    end

    curDataSet.add(wekaInstance);
    wekaPredictedClassLabel = wekaModelClassifier.classifyInstance(curDataSet.instance(0)) + 1;
    predictedLabel = wekaModelClassNameList{wekaPredictedClassLabel};
    classPredictionProbabilities = wekaModelClassifier.distributionForInstance( curDataSet.instance(0) );

end