function [ cc ] = ccAppend( cc, pixels )
%ccAppend Append a new connected component to the list

if(length(pixels) ~= numel(pixels))
    pixels = sub2ind(cc.ImageSize,pixels(:,1),pixels(:,2));
end

cc.NumObjects = cc.NumObjects + 1;
cc.PixelIdxList{cc.NumObjects} = pixels(:);


end

