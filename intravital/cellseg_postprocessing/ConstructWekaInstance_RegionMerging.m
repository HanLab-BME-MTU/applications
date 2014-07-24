function [wekaInstance, wekaAttributeList] = ConstructWekaInstance_RegionMerging(featureVec, featureNameList) 

    classNameList = {'Dont-Merge', 'Merge'};
    
    wekaClassNameList = weka.core.FastVector();
    for i = 1:numel(classNameList)
        wekaClassNameList.addElement( java.lang.String(classNameList{i}) );
    end
    
    wekaAttributeList = weka.core.FastVector();                                
    wekaFeatureVec = zeros(1, numel(featureNameList));

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
