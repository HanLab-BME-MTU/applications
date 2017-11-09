function bounded = ccBound(cc)
    bounded = cc;
    for i=1:cc.NumObjects
        cc.PixelIdxList{i} = i;
    end
end
