function saveFigure3(inputDirectory, outputDirectory)

Z = zeros(400, 400, 'uint8');

paths = {'TM2_Actin/cell5_Q4', 'TM4_Actin/14_August_2009/cell2', 'TM5NM1_Actin/26June2009/Cell6_Q4'};
actinPaths = {'Actin', 'actin', 'Actin'};
tmPaths = {'TM2', 'TM4', 'TM5NM1'};
imageFileName = {'default01.tif', 'crop_default002.tif', 'default080.tif'};
pos = [55, 149; 97, 115; 62, 23];

for i = 1:3
    %
    % Panel A
    %
    
    % Actin
    path = [inputDirectory filesep paths{i} filesep actinPaths{i}];
    filename = [path filesep 'crop' filesep imageFileName{i}];
    I1 = imread(filename);
    % Crop
    I1 = I1(pos(i,1):pos(i,1)+400-1,pos(i,2):pos(i,2)+400-1);
    % Convert to 8bit channel
    iMin = min(I1(:));
    iMax = max(I1(:));
    I1 = uint8((255 / double(iMax - iMin)) * (double(I1) - iMin));
    C = cat(3, Z, I1, Z);
    %imwrite(C, [outputDirectory filesep 'fig3_A' num2str(i) '1.tif'], 'Compression', 'none');
    
    % TM
    path = [inputDirectory filesep paths{i} filesep tmPaths{i}];
    filename = [path filesep 'crop' filesep imageFileName{i}];
    I2 = imread(filename);
    % Crop
    I2 = I2(pos(i,1):pos(i,1)+400-1,pos(i,2):pos(i,2)+400-1);
    % Convert to 8bit channel
    iMin = min(I2(:));
    iMax = max(I2(:));
    I2 = uint8((255 / double(iMax - iMin)) * (double(I2) - iMin));
    C = cat(3, I2, Z, Z);
    %imwrite(C, [outputDirectory filesep 'fig3_A' num2str(i) '2.tif'], 'Compression', 'none');
    
    % Merge
    C = cat(3, I1, I2, Z);
    %imwrite(C, [outputDirectory filesep 'fig3_A' num2str(i) '3.tif'], 'Compression', 'none');
    
    %
    % Panel B
    %
    hFig = figure;
    imshow(C);
end

%
% Panel B
%