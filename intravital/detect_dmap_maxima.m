function [ maxpts ] = detect_dmap_maxima( imdist )

    imsize = size( imdist );
    
    maxpts = [];
    for i = 2:imsize(1)-1
        for j = 2:imsize(2)-1
            
            neighVals = [ imdist(i-1,j), imdist(i+1,j), imdist(i,j-1), imdist(i,j+1) ];
            
            if imdist(i,j) > max(neighVals)
                maxpts(end+1,:) = [j,i];
            end
            
        end
    end

end