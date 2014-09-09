function dilatedCC = ccDilate(cc,se)
    import connectedComponents.*;
    if(isa(se,'strel'))
        se = se.getnhood;
    end
    
    dilatedCC = ccPad(cc,size(se));
    offsets = mat2offset(dilatedCC.ImageSize,se);


    for i=1:cc.NumObjects
        dilatedCC.PixelIdxList{i} = bsxfun(@plus,dilatedCC.PixelIdxList{i},offsets);
        dilatedCC.PixelIdxList{i} = unique(dilatedCC.PixelIdxList{i});
        dilatedCC.PixelIdxList{i} = dilatedCC.PixelIdxList{i}(:);
    end
    dilatedCC = ccUnpad(dilatedCC);
end
