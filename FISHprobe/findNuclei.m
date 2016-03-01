function nuclei = findNuclei(handles)
%findNuclei Crop out all nuclei with three different channels
%   Detailed explanation goes here

% To be optimized
% handles is a parameter from InvivoCytometer_2.0_source_code/code_package/CellSegmentationQualityAnnotator.m

% 02/2016 Ning

dapiCha = input('Enter the DAPI channel number > ');
greenCha = input('Enter the Green channel number > ');
redCha = input('Enter the Red channel number > ');

mask = handles.data.imLabelCellSeg;
numNuclei = max(mask(:));

% Try to pre allocate the nuclei structure size

for num = 1:numNuclei
    nucMask = zeros(size(mask));
    nucMask(mask==num) = 1;
    
    
    % Modify the dilation factors as is
    % se=strel('ball',10,3);
    [a,b,c] = ndgrid(-10:10,-10:10,-3:3);
    se = strel(sqrt(a.^2+b.^2+c.^2));
    nucMask2 = imdilate(nucMask,se);

    nucDapi = handles.data.imageData{dapiCha};    
    nucDapi(nucMask2==0) = 0;
    
    nucGreen = handles.data.imageData{greenCha};
    nucRed = handles.data.imageData{redCha};
    
    lowerX = 1;
    lowerY = 1;
    lowerZ = 1;
    
    for i = 1:size(mask,1)
        if min(min(nucDapi(i,:,:)==0))==0
            if lowerX == 1
                lowerX = i;
            end
            higherX = i;
        end
    end
    
    for j = 1:size(mask,2)
        if min(min(nucDapi(:,j,:)==0))==0
            if lowerY == 1
                lowerY = j;
            end
            higherY = j;
        end
    end

    for k = 1:size(mask,3)
        if min(min(nucDapi(:,:,k)==0))==0
            if lowerZ == 1
                lowerZ = k;
            end
            higherZ = k;
        end
    end
    
    
    % How to predefine structure size???
    nuclei(num).dapi = nucDapi(lowerX:higherX, lowerY:higherY, lowerZ:higherZ);
    nuclei(num).green = nucGreen(lowerX:higherX, lowerY:higherY, lowerZ:higherZ);
    nuclei(num).red = nucRed(lowerX:higherX, lowerY:higherY, lowerZ:higherZ);
end


end

