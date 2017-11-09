function imageWithNumber=addNumbersImage(imagen, number, numberSize, position)

% adds to an image a legend with a number embedded into it
% Example
% addNumbersImage(imagen, 146, [25, 75], [20, 10]) will embed the number 146
% into the image with a size of [25 75] with top right corner [20, 10]

% Author: Santiago Costantino 
% santiago.costantino@umontreal.ca

%%Check the numbers image fits in the image, if not, shits it until it does

warning('off')

if numberSize(1)+position(1)-1>size(imagen, 1)
    position(1)=position(1)-numberSize(1);
end

if numberSize(2)+position(2)-1>size(imagen, 2)
    position(2)=position(2)-numberSize(2);
end

%% Adds the number
if size(imagen, 3)<3
    numberImage=zeros(size(imagen));
    numberImage(position(1):numberSize(1)+position(1)-1, ...
        position(2):numberSize(2)+position(2)-1)=makeCharaterImage(number, numberSize)*255;
    imagen(numberImage~=0)=numberImage(numberImage~=0);
    imageWithNumber=imagen;
else
    for itColor=1:3
        thisColor=imagen(:,:,itColor);
        numberImage=zeros(size(thisColor));
        numberImage(position(1):numberSize(1)+position(1)-1, ...
            position(2):numberSize(2)+position(2)-1)=makeCharaterImage(number, numberSize)*255;
        thisColor(numberImage~=0)=numberImage(numberImage~=0);
        imageWithNumber(:,:,itColor)=thisColor;
    end
end

warning('on')