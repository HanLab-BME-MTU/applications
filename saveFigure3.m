function saveFigure3(inputDirectory, outputDirectory)

paths = {'TM2_Actin/cell5_Q4', 'TM4_Actin/14_August_2009/cell2', 'TM5NM1_Actin/26June2009/Cell6_Q4'};
actinPaths = {'Actin', 'actin', 'Actin'};
tmPaths = {'TM2', 'TM4', 'TM5NM1'};
imageFileName = {'default01.tif', 'crop_default002.tif', 'default080.tif'};
locMaxFileName = {'locMax01.mat', 'locMax002.mat', 'locMax080.mat'};
w = 400;
pos = [55, 149; 97, 115; 62, 23];
wInset = 50;
posInset = [];

Z = zeros(w, w, 'uint8');

for i = 1:3
    %
    % Panel A
    %
    
    % Actin
    path = [inputDirectory filesep paths{i} filesep actinPaths{i}];
    filename = [path filesep 'crop' filesep imageFileName{i}];
    I1 = imread(filename);
    % Crop
    I1 = I1(pos(i,1):pos(i,1)+w-1,pos(i,2):pos(i,2)+w-1);
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
    I2 = I2(pos(i,1):pos(i,1)+w-1,pos(i,2):pos(i,2)+w-1);
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
    path = [inputDirectory filesep paths{i} filesep actinPaths{i}];
    filename = [path filesep 'edge' filesep 'cell_mask' filesep 'mask_' imageFileName{i}];
    Bw = imread(filename);
    % Crop
    Bw = Bw(pos(i,1):pos(i,1)+w-1,pos(i,2):pos(i,2)+w-1);
    D = double(bwdist(1 - Bw)) * 67;
    %A = NaN(size(Bw));
    %A(D < 5000 & Bw == 1) = 1;
    figure, imshow(C);
    % Add Actin speckles
    path = [inputDirectory filesep paths{i} filesep actinPaths{i}];
    filename = [path filesep 'tack' filesep 'locMax' filesep locMaxFileName{i}];
    load(filename);
    % Crop
    locMax = locMax(pos(i,1):pos(i,1)+w-1,pos(i,2):pos(i,2)+w-1);
    idx = find(locMax ~= 0 & D < 5000);
    [y x] = ind2sub(size(locMax), idx);
    line(x, y,'LineStyle', 'none', 'Marker', '.', 'Color', 'r','MarkerSize',6);
    % Add TM speckles
    path = [inputDirectory filesep paths{i} filesep tmPaths{i}];
    filename = [path filesep 'tack' filesep 'locMax' filesep locMaxFileName{i}];
    load(filename);
    % Crop
    locMax = locMax(pos(i,1):pos(i,1)+w-1,pos(i,2):pos(i,2)+w-1);
    idx = find(locMax ~= 0 & D < 5000);
    [y x] = ind2sub(size(locMax), idx);
    line(x, y,'LineStyle', 'none', 'Marker', '.', 'Color', 'g','MarkerSize',6);
    c = contourc(double(D), [0; 5000]);
    n = c(2, 1);
    line(c(1, 2:n+2), c(2, 2:n+2), 'Color', 'w', 'Linewidth', .75);
    line(c(1, n+3:end), c(2, n+3:end), 'Color', 'w', 'Linewidth', .75);
    print(gcf, '-depsc' , '-painters', [outputDirectory filesep 'fig3_B' num2str(i) '.eps']);
end