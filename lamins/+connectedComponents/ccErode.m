function cc = ccErode(cc,se)
    L = labelmatrix(cc);
    bg = imdilate(~L,se);
    cc.PixelIdxList = cellfun(@(x) x(~bg(x)),cc.PixelIdxList,'UniformOutput',false);
end
