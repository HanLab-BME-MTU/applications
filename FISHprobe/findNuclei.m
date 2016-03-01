function nuclei = findNuclei(imageData, mask)
%findNuclei Crop out all nuclei with three different channels
%   Detailed explanation goes here

% To be optimized
% imageData = handles.data.imageData; (Multi-channel 3D stack)
% mask = handles.data.imLabelCellSeg; (3D stack for Dapi channel)
% handles is a parameter from InvivoCytometer_2.0_source_code/code_package/CellSegmentationQualityAnnotator.m

% 03/2016 Ning

dapiCha = input('Enter the DAPI channel number > ');
greenCha = input('Enter the Green channel number > ');
redCha = input('Enter the Red channel number > ');

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

    nucDapi = imageData{dapiCha};    
    nucDapi(nucMask2==0) = 0;
    
    nucGreen = imageData{greenCha};
    nucRed = imageData{redCha};
    
    xStart = 1;
    yStart = 1;
    zStart = 1;
    
    for i = 1:size(mask,1)
        if min(min(nucDapi(i,:,:)==0))==0
            if xStart == 1
                xStart = i;
            end
            xEnd = i;
        end
    end
    
    for j = 1:size(mask,2)
        if min(min(nucDapi(:,j,:)==0))==0
            if yStart == 1
                yStart = j;
            end
            yEnd = j;
        end
    end

    for k = 1:size(mask,3)
        if min(min(nucDapi(:,:,k)==0))==0
            if zStart == 1
                zStart = k;
            end
            zEnd = k;
        end
    end
    
    
    % How to predefine structure size???
    nuclei(num).dapi = nucDapi(xStart:xEnd, yStart:yEnd, zStart:zEnd);
    nuclei(num).green = nucGreen(xStart:xEnd, yStart:yEnd, zStart:zEnd);
    nuclei(num).red = nucRed(xStart:xEnd, yStart:yEnd, zStart:zEnd);
    nuclei(num).xRange = [xStart, xEnd];
    nuclei(num).yRange = [yStart, yEnd];
    nuclei(num).zRange = [zStart, zEnd];
end


end

