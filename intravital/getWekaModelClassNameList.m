function [wekaModelClassNameList] = getWekaModelClassNameList( wekaModelAttributeHeader )

   wekaModelClassNameList = cell(1, wekaModelAttributeHeader.attribute(0).numValues);
   for i = 1:numel(wekaModelClassNameList)
       wekaModelClassNameList{i} = char( wekaModelAttributeHeader.attribute(0).value(i-1) );
   end

end