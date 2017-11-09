function cc = ccBinaryOp(func,A,B)
    assert(A.NumObjects == B.NumObjects);
    assert(all(A.ImageSize == B.ImageSize));
    cc = A;
    cc.PixelIdxList = cellfun(func,A.PixelIdxList,B.PixelIdxList,'UniformOutput',false);
end
