function [flagEqual, strErrorMessage] = wekaCheckIfHeadersAreEqual(wekaModelHeader, wekaDataset)

    % check if class indices are equal
    if wekaModelHeader.classIndex() ~= wekaDataset.classIndex()       
        flagEqual = false;
        strErrorMessage = sprintf( 'Class index diffex: %d ~= %d', ...
                                    wekaModelHeader.classIndex(), ...
                                    wekaDataset.classIndex() );
        return;
    end
    
    % check if size of the attribute vectors are equal
    if wekaModelHeader.numAttributes() ~= wekaDataset.numAttributes()       
        flagEqual = false;
        strErrorMessage = sprintf( 'Different number of attributes: %d ~= %d', ...
                                    wekaModelHeader.numAttributes(), ... 
                                    wekaDataset.numAttributes() );
        return;
    end
    
    % check if all the attribute names are equal
    for i = 1:wekaModelHeader.numAttributes()
        
        modelAttribute = wekaModelHeader.attribute(i-1);
        datasetAttribute = wekaDataset.attribute(i-1);
        
        % check if names are equal
        if ~strcmp( char(modelAttribute.name()), ...
                    char(datasetAttribute.name()) )
                
            flagEqual = false;
            strErrorMessage = sprintf( 'Attribute names differ at position %d: %s ~= %s', ...
                                        i, char(modelAttribute.name()), ...
                                        char(datasetAttribute.name()) );
            return;
        end
        
        % check types
        if modelAttribute.type() ~= datasetAttribute.type()
            
            flagEqual = false;
            
            type1 = modelAttribute.type();
            type2 = datasetAttribute.type();
            strType1 = char(weka.core.Attribute.typeToString(type1));
            strType2 = char(weka.core.Attribute.typeToString(type2));
            
            strErrorMessage = sprintf( 'Types of attribute %s at position %d differ: %s ~= %s', ...
                                        char(modelAttribute.name()), i, ...
                                        strType1, strType2 );
            return;
        end

        % if nominal check if the values are same
        if modelAttribute.isNominal() && datasetAttribute.isNominal()
            
            if modelAttribute.numValues() ~= datasetAttribute.numValues()
               flagEqual = false;
               
                strErrorMessage = sprintf( 'Number of values of nominal attribute %s at position %d differ: %d ~= %d', ...
                                            char(modelAttribute.name()), i, ...
                                            modelAttribute.numValues(), ...
                                            datasetAttribute.numValues() );
               
               return;
            end

            for vid = 1:modelAttribute.numValues()
                
                if ~strcmp( char(modelAttribute.value(vid-1)), ...
                            char(datasetAttribute.value(vid-1)) )

                    flagEqual = false;
                        
                    modelAttributeValues = [];
                    datasetAttributeValues = [];
                    
                    for vid2 = 1:modelAttribute.numValues()
                        
                        curPrefix = '';
                        if vid2 > 1
                            curPrefix = ', '; 
                        end
                        
                        modelAttributeValues = [modelAttributeValues, curPrefix, char(modelAttribute.value(vid2-1))];
                        datasetAttributeValues = [datasetAttributeValues, curPrefix, char(datasetAttribute.value(vid2-1))];
                        
                    end
                    
                    strErrorMessage = sprintf( 'Values of nominal attribute %s at position %d differ: [%s] ~= [%s]', ...
                                                char(modelAttribute.name()), i, ...
                                                modelAttributeValues, datasetAttributeValues);
                    return;                
                    
                end
                
            end
            
        end
        
    end

    flagEqual = true;
    strErrorMessage = [];
    
end