function nucleiStruc = spotsPair(nucleiStruc, dataProperties)
%spotsPair Finds and pairs the nearest spot in the other channel;
%calculates the mutual distance
%   Detailed explanation goes here

% Get pixel size in um.
pixelSizeXY = dataProperties.PIXELSIZE_XY;
pixelSizeZ = dataProperties.PIXELSIZE_Z;

for nucNum = 1:numel(nucleiStruc)    
    redSpotsPool = nucleiStruc(nucNum).redSpot;
    greenSpotsPool = nucleiStruc(nucNum).greenSpot;
    redSpotPair = {};
    greenSpotPair = {};
    distance = [];
    pairCount = 0;
    
    for i = 1:numel(redSpotsPool)
        redSpot = redSpotsPool(i).cord;
        pairFound = 0;
        distLimit = 3;
        for j = 1:numel(greenSpotsPool)
            greenSpot = greenSpotsPool(j).cord;
            dist = ( ((greenSpot(1) - redSpot(1)) * pixelSizeXY)^2 ...
                + ((greenSpot(2) - redSpot(2)) * pixelSizeXY)^2 ...
                + ((greenSpot(3) - redSpot(3)) * pixelSizeZ)^2 ) ^ (1/2);
            
            % Remove spots with pair distance larger than 3um;
            % Choose nearest spot pair if two or more are within 3um
            % May use smarter cluster method
            
            if dist > distLimit
                continue;
            else
                pairFound = 1;
                distLimit = dist;
                greenSpotRecord = greenSpot;
            end            
        end
        
        if pairFound
            % Check if one green dot is paired with multiple red dots
            pairDuplicate = 0;
            for spotNum = 1:pairCount
                if any(greenSpotRecord == greenSpotPair(spotNum).cord)
                    pairDuplicate = 1;
                    if distLimit < distance(spotNum)
                        redSpotPair(spotNum).cord = redSpot;
                    end
                end
            end
            
            if ~pairDuplicate
                pairCount = pairCount +1 ;
                redSpotPair(pairCount).cord = redSpot;
                greenSpotPair(pairCount).cord = greenSpotRecord;
                distance(pairCount) = distLimit;
            end
        end
        
    end
    nucleiStruc(nucNum).redSpotPair = redSpotPair;
    nucleiStruc(nucNum).greenSpotPair = greenSpotPair;
    nucleiStruc(nucNum).distance = distance;
    
end

end

