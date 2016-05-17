function nucleiStruc = findNuclei(imageData, mask, dataProperties)
%findNuclei Crop out all nuclei with three different channels
%   Detailed explanation goes here

% 04/2016 Ning Zhang

numNuclei = max(mask(:));

% Preallocate nuclei structure
nucleiStruc(numNuclei).dapi = [];
for chaNum = 1:numel(dataProperties.channel)
    chaName = dataProperties.channel(chaNum).name;
    switch chaName
        case 'dapi'
            continue;
        case 'green' 
            nucleiStruc(numNuclei).green = [];
        case 'red'
            nucleiStruc(numNuclei).red = [];
        otherwise
            error('Unknown channels detected')
    end
end
nucleiStruc(numNuclei).xRange = [];
nucleiStruc(numNuclei).yRange = [];
nucleiStruc(numNuclei).zRange = [];

% Crop single nucleus one by one
for num = 1:numNuclei
    nucMask = zeros(size(mask));
    nucMask(mask==num) = 1;

    % Modify the dilation factors as is
    % se=strel('ball',10,3);
%     [a,b,c] = ndgrid(-10:10,-10:10,-3:3);
%     se = strel(sqrt(a.^2+b.^2+c.^2));
%     nucMask = imdilate(nucMask,se);
    % Make imdilate faster!!!

    nucDapi = imageData.dapi;    
    nucDapi(nucMask==0) = 0;

    % Inverted x, y for matrix and cartesian coordinates
    yStart = 1;
    xStart = 1;
    zStart = 1;

    for i = 1:size(mask,1)
        % when all layer values equal to zero is not true
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
    
    for chaNum = 1:numel(dataProperties.channel)
        chaName = dataProperties.channel(chaNum).name;
        switch chaName
            case 'dapi'
                continue;
            case 'green' 
                nucGreen = imageData.green;
                nucGreen(nucMask==0) = 0;
                nucleiStruc(num).green = nucGreen(yStart:yEnd, xStart:xEnd, zStart:zEnd);
                clear nucGreen
            case 'red'
                nucRed = imageData.red;
                nucRed(nucMask==0) = 0;
                nucleiStruc(num).red = nucRed(yStart:yEnd, xStart:xEnd, zStart:zEnd);
                clear nucRed
            otherwise
                error('Unknown channels detected')
        end
    end
    nucleiStruc(num).xRange = [xStart, xEnd];
    nucleiStruc(num).yRange = [yStart, yEnd];
    nucleiStruc(num).zRange = [zStart, zEnd];
    
    clear nucMask nucDapi
end

end

