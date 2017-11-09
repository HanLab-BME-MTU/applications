function [ angularLocalMaxima ] = findRobustPeaksOnCircle( angularResponse )
%findRobustPeaksOnCircle Finds peaks on a circle where 

    angularLocalMaxima = angularResponse > angularResponse(:,:,[end 1:end-1]) ...
                       & angularResponse > angularResponse(:,:,[2:end 1]) ...
                       & angularResponse(:,:,[2:end 1]) > angularResponse(:,:,[3:end 1 2]) ...
                       & angularResponse(:,:,[end 1:end-1]) > angularResponse(:,:,[end-1 end 1:end-2]);


end

