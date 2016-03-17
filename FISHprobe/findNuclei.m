function nucleiStruc = findNuclei(imageData, mask)
%findNuclei Crop out all nuclei with three different channels
%   Detailed explanation goes here

% 03/2016 Ning

numNuclei = max(mask(:));

% Preallocate nuclei structure
nucleiStruc(numNuclei).dapi = [];
nucleiStruc(numNuclei).green = [];
nucleiStruc(numNuclei).red = [];
nucleiStruc(numNuclei).xRange = [];
nucleiStruc(numNuclei).yRange = [];
nucleiStruc(numNuclei).zRange = [];


% Try to pre allocate the nuclei structure size

for num = 1:numNuclei
    nucMask = zeros(size(mask));
    nucMask(mask==num) = 1;
    
    
    % Modify the dilation factors as is
    % se=strel('ball',10,3);
%     [a,b,c] = ndgrid(-10:10,-10:10,-3:3);
%     se = strel(sqrt(a.^2+b.^2+c.^2));
%     nucMask2 = imdilate(nucMask,se);
    % Make imdilate faster!!!
    
    nucMask2 = nucMask;

    nucDapi = imageData.dapi;    
    nucDapi(nucMask2==0) = 0;
    
    nucGreen = imageData.green;
    nucRed = imageData.red;
    
    % Inverted x, y for matrix and coordinates
    
    yStart = 1;
    xStart = 1;
    zStart = 1;
    
    for i = 1:size(mask,1)
        if min(min(nucDapi(i,:,:)==0))==0
            if yStart == 1
                yStart = i;
            end
            yEnd = i;
        end
    end
    
    for j = 1:size(mask,2)
        if min(min(nucDapi(:,j,:)==0))==0
            if xStart == 1
                xStart = j;
            end
            xEnd = j;
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
    
    nucleiStruc(num).dapi = nucDapi(yStart:yEnd, xStart:xEnd, zStart:zEnd);
    nucleiStruc(num).green = nucGreen(yStart:yEnd, xStart:xEnd, zStart:zEnd);
    nucleiStruc(num).red = nucRed(yStart:yEnd, xStart:xEnd, zStart:zEnd);
    nucleiStruc(num).xRange = [xStart, xEnd];
    nucleiStruc(num).yRange = [yStart, yEnd];
    nucleiStruc(num).zRange = [zStart, zEnd];
end


end

