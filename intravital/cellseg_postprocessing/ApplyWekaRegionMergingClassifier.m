function [predictedClassLabel, predictionProbability] = ApplyWekaRegionMergingClassifier(wekaModel, featureVec, featureNameList)

    import weka.core.*;
    
    wekaModelClassifier = wekaModel(1);
    wekaModelAttributeHeader = wekaModel(2);
    
    wekaModelClassNameList = cell(1, wekaModelAttributeHeader.attribute(0).numValues);
    for i = 1:numel(wekaModelClassNameList)
        wekaModelClassNameList{i} = char( wekaModelAttributeHeader.attribute(0).value(i-1) );
    end
    
    [wekaInstance, wekaAttributeList] = ConstructWekaInstance_RegionMerging(featureVec, featureNameList, wekaModelClassNameList);

    curDataSet = weka.core.Instances('TestInstance', wekaAttributeList, 1 );
    curDataSet.add(wekaInstance);
    curDataSet.setClassIndex(0);

    if ~wekaModelAttributeHeader.equalHeaders(curDataSet)
        error( 'ERROR: mismatch between computed features and classifier features' );
    end

    wekaPredictedClassLabel = wekaModelClassifier.classifyInstance(curDataSet.instance(0)) + 1;
    wekaPredictionProbabilityDist = wekaModelClassifier.distributionForInstance( curDataSet.instance(0) );
    
    switch wekaModelClassNameList{wekaPredictedClassLabel}
    
        case 'Merge'
            
            predictedClassLabel = 1;
            
        case 'Dont-Merge'
        
            predictedClassLabel = 0;
            
    end    
    
    predictionProbability = wekaPredictionProbabilityDist(wekaPredictedClassLabel);
    
end

function [wekaInstance, wekaAttributeList] = ConstructWekaInstance_RegionMerging(featureVec, featureNameList, classNameList) 

    import weka.core.*;

    wekaAttributeList = weka.core.FastVector();                                
    wekaFeatureVec = zeros(1, numel(featureNameList));

    wekaClassNameList = weka.core.FastVector();
    for i = 1:numel(classNameList)
        wekaClassNameList.addElement( java.lang.String(classNameList{i}) );
    end
    
    classAttributeName = java.lang.String('ClassLabel');
    classAttribute = weka.core.Attribute(classAttributeName, wekaClassNameList);
    classAttributeVal = 0; % dummy label
    wekaAttributeList.addElement(classAttribute);
    wekaFeatureVec(1) = classAttributeVal;
    
    for attid = 1:numel(featureNameList)                                                                        

        curAttributeName = java.lang.String( featureNameList{attid} );

        curAttribute = weka.core.Attribute(curAttributeName);
        curAttributeVal = featureVec{attid};

        wekaAttributeList.addElement(curAttribute);
        wekaFeatureVec(1+attid) = curAttributeVal;

    end

    wekaInstance = weka.core.DenseInstance(1.0, wekaFeatureVec);

end
