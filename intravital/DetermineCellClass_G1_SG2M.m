function [cellFeatureClass, classNameList] = DetermineCellClass_G1_SG2M( inCellPatternType )

    cellFeatureClass = [];
    classNameList = { 'G1', 'S-G2-M' };

    switch inCellPatternType

        case { 'Mono_Pre_G1', 'Multi_Pre_G1', 'Mono_G1', 'Multi_G1', 'G1' } 

            cellFeatureClass = 'G1';

        case { 'Mono_S', 'Multi_S', 'S' } 

            cellFeatureClass = 'S-G2-M';

        case { 'Mono_G2', 'Multi_G2', 'G2' }

            cellFeatureClass = 'S-G2-M';

        case { 'Mono_Prophase', 'Multi_Prophase', 'Mono_Anaphase', 'Multi_Anaphase', 'Mono_Metaphase', 'Multi_Metaphase', 'M' }     

            cellFeatureClass = 'S-G2-M';

    end   
    
end