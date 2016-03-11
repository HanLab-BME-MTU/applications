function position=getCenterOfMass(imagen, binaryMask)


% Author: Santiago Costantino 
% santiago.costantino@umontreal.ca


% outputs the center of mass of an image, considering only the pixels in
% the mask

maskedImage=double(imagen).*double(binaryMask);
intensitySum=sum(sum(maskedImage));

position=[0, 0];
for itRow=1:size(imagen, 1)
    for itColumn=1:size(imagen, 2)
        position=position+[maskedImage(itRow, itColumn)*itRow, maskedImage(itRow, itColumn)*itColumn];
    end
end
position=position./intensitySum;
%% Change to x,y coordinates
position=[position(2) position(1)];