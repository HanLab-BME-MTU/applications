function images = getImages(obj)
    import lamins.classes.LaminsImage;
    for c = 1 : obj.reader.getSizeC
        for t = 1 : obj.reader.getSizeT
            for z = 1 : obj.reader.getSizeZ
                images(c,t,z) = LaminsImage(obj,c,t,z);
            end
        end
    end
end
