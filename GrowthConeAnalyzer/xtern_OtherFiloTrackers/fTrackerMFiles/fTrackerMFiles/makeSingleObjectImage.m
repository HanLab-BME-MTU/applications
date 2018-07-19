function singleObject=makeSingleObjectImage(multipleObjects)

% it deletes from an image all the objects but the bigest.

% Author: Santiago Costantino 
% santiago.costantino@umontreal.ca


[multipleObjectsLabeled, labels]=bwlabel(multipleObjects);
objectAreas=regionprops(multipleObjectsLabeled, 'MajorAxisLength', 'MinorAxisLength');

maximo=1;
if max(labels>1)
    for it=1:max(labels);
        if maximo<(objectAreas(it,1).MajorAxisLength)*(objectAreas(it,1).MinorAxisLength)
            maximo=(objectAreas(it,1).MajorAxisLength)*(objectAreas(it,1).MinorAxisLength);
            label=it;
        end
    end
else
    label=1;
end
        
singleObject=zeros(size(multipleObjects));
singleObject(multipleObjectsLabeled==label)=1;