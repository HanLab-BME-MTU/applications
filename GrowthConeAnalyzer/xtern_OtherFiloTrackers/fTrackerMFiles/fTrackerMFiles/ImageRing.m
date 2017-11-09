function mat = ImageRing(sizeIm,radius,x,y,value,width)

% Author: Antoine Godin
% godin.antoine@sympatico.ca

sizeX = sizeIm(1);
switch length(sizeIm)
    case 1 
        sizeY = sizeX;
    case 2
        sizeY = sizeIm(2);
end

if x < radius+1
    xmin = 1;
    xmax = 2*radius;
else
    if x > (sizeX-radius+1)
        xmin = sizeX - 2*radius+1;
        xmax = sizeX;
    else
        xmin = (x-(radius+1)+1);
        xmax = (x+(radius+1));
    end
end
if y < radius+1
    ymin = 1;
    ymax = 2*radius;
else
    if y > (sizeY-radius+1)
        ymin = sizeY - 2*radius+1;
        ymax = sizeY;
    else
        ymin = (y-(radius+1)+1);
        ymax = (y+(radius+1));
    end
end


if value == 1
    mat = zeros(sizeIm);
    for itx = xmin:min(xmax,sizeX)
        for ity = ymin:min(ymax,sizeY)
            if ((itx-x)^2 + (ity-y)^2 ) <= (radius+width/2)^2 & ...
                    ((itx-x)^2 + (ity-y)^2 ) >= (radius-width/2)^2
                mat(itx,ity) = 1;
            end
        end
    end
else
    mat = ones(sizeIm);
    for itx = xmin:min(xmax,sizeX)
        for ity = ymin:min(ymax,sizeY)
            if ((itx-x)^2 + (ity-y)^2 ) <= (radius+width/2)^2 & ...
                    ((itx-x)^2 + (ity-y)^2 ) >= (radius-width/2)^2
               mat(itx,ity) = 0;
            end
        end
    end
end