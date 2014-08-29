function mask = bigmask(L)
%     L = Lamins(1);
    Lcell = squeeze(L.cellReader.toCell);
    mask = cell(size(Lcell));
    for c = 1:size(Lcell,1);
        parfor z = 1:size(Lcell,2);
            mask{c,z} = LaminsImage(Lcell{c,z}).mask;
        end;
    end;
end