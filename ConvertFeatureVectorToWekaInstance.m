function [wekaInstance, wekaAttributeList, classNameList] = ConvertFeatureVectorToWekaInstance(wekaModel, featureVec, featureNameList)

    wekaModelAttributeHeader = wekaModel(2);
    wekaModelClassAttribute = wekaModelAttributeHeader.attribute(0);

    if ~(wekaModelClassAttribute.isNominal() && ...
      strcmp( char(wekaModelClassAttribute.name()), 'ClassLabel' ))
      error('The first feature must be ClassLabel'); 
    end

    wekaClassNameList = weka.core.FastVector();
    classNameList = cell(1, wekaModelClassAttribute.numValues);
    for i = 1:numel(classNameList)
       classNameList{i} = char( wekaModelClassAttribute.value(i-1) );
       wekaClassNameList.addElement( java.lang.String(classNameList{i}) );
    end

    wekaFeatureVec = zeros(1, numel(featureNameList));
    wekaAttributeList = weka.core.FastVector();

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
