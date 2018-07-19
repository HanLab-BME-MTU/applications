function [cellFeatureClass, classNameList] = DetermineCellClass_Interphase_Mitotic( inCellPatternType )

    cellFeatureClass = [];
    classNameList = { 'Interphase', 'Mitotic' };
    
    switch inCellPatternType

        case { 'Mono_Pre_G1', 'Multi_Pre_G1', 'Mono_G1', 'Multi_G1', 'G1' } 

            cellFeatureClass = 'Interphase';

        case { 'Mono_S', 'Multi_S', 'S' } 

            cellFeatureClass = 'Interphase';

        case { 'Mono_G2', 'Multi_G2', 'G2' }

            cellFeatureClass = 'Interphase';

        case { 'Mono_Prophase', 'Multi_Prophase', 'Mono_Anaphase', 'Multi_Anaphase', 'Mono_Metaphase', 'Multi_Metaphase', 'M' }     

            cellFeatureClass = 'Mitotic';

    end    

end