function [cellFeatureClass, classNameList] = DetermineCellClass_S_G2_M( inCellPatternType )

    cellFeatureClass = [];
    classNameList = { 'S', 'G2', 'M' };
    
    switch inCellPatternType

        case { 'Mono_S', 'Multi_S', 'S' } 

            cellFeatureClass = 'S';

        case { 'Mono_G2', 'Multi_G2', 'G2' }

            cellFeatureClass = 'G2';

        case { 'Mono_Prophase', 'Multi_Prophase', 'M' }     

            cellFeatureClass = 'M';

    end    

end