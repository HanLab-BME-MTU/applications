function [cellFeatureClass, classNameList] = DetermineCellClass_G2_M( inCellPatternType )

    cellFeatureClass = [];
    classNameList = { 'G2', 'M' };
    
    switch inCellPatternType

        case { 'Mono_G2', 'Multi_G2', 'G2' }

            cellFeatureClass = 'G2';

        case { 'Mono_Prophase', 'Multi_Prophase', 'Mono_Anaphase', 'Multi_Anaphase', 'Mono_Metaphase', 'Multi_Metaphase', 'M' }     

            cellFeatureClass = 'M';

    end    

end