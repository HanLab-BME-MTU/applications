function saveFigure3(inputDirectory, outputDirectory)

paths = {'TM2_Actin/cell5_Q4', 'TM4_Actin/14_August_2009/cell2', 'TM5NM1_Actin/26June2009/Cell6_Q4'};
actinPaths = {'Actin', 'actin', 'Actin'};
tmPaths = {'TM2', 'TM4', 'TM5NM1'};
imageFileName = {'default01.tif', 'crop_default002.tif', 'default080.tif'};
locMaxFileName = {'locMax01.mat', 'locMax002.mat', 'locMax080.mat'};
w = 400;
pos = [55, 149; 97, 115; 62, 23];
wInset = 50;
posInset = [200,285; 246, 204; 372, 279];

Z1 = zeros(w, w, 'uint8');
Z2 = zeros(wInset, wInset, 'uint8');

for i = 1:3
    %
    % Panel A
    %
    
    % Actin
    path = [inputDirectory filesep paths{i} filesep actinPaths{i}];
    filename = [path filesep 'crop' filesep imageFileName{i}];
    IActin = imread(filename);
    % Convert to 8bit channel
    iMin = min(IActin(:));
    iMax = max(IActin(:));
    IActin = uint8((255 / double(iMax - iMin)) * (double(IActin) - iMin));
    % Crop
    IActin1 = IActin(pos(i,1):pos(i,1)+w-1,pos(i,2):pos(i,2)+w-1);
    C1 = cat(3, IActin1, Z1, Z1);
    %imwrite(C1, [outputDirectory filesep 'fig3_A' num2str(i) '1.tif'], 'Compression', 'none');
    
    % TM
    path = [inputDirectory filesep paths{i} filesep tmPaths{i}];
    filename = [path filesep 'crop' filesep imageFileName{i}];
    ITM = imread(filename);
    % Convert to 8bit channel
    iMin = min(ITM(:));
    iMax = max(ITM(:));
    ITM = uint8((255 / double(iMax - iMin)) * (double(ITM) - iMin));
    % Crop
    ITM1 = ITM(pos(i,1):pos(i,1)+w-1,pos(i,2):pos(i,2)+w-1);
    C1 = cat(3, Z1, ITM1, Z1);
    %imwrite(C1, [outputDirectory filesep 'fig3_A' num2str(i) '2.tif'], 'Compression', 'none');
    
    % Merge
    C1 = cat(3, IActin1, ITM1, Z1);
    %imwrite(C1, [outputDirectory filesep 'fig3_A' num2str(i) '3.tif'], 'Compression', 'none');
    
    %
    % Panel B
    %
    path = [inputDirectory filesep paths{i} filesep actinPaths{i}];
    filename = [path filesep 'edge' filesep 'cell_mask' filesep 'mask_' imageFileName{i}];
    Bw = imread(filename);
    % Compute distance transform
    D = double(bwdist(1 - Bw)) * 67;
    % Load actin speckles
    path = [inputDirectory filesep paths{i} filesep actinPaths{i}];
    filename = [path filesep 'tack' filesep 'locMax' filesep locMaxFileName{i}];
    load(filename);
    locActin = locMax;
    % Load TM speckles
    path = [inputDirectory filesep paths{i} filesep tmPaths{i}];
    filename = [path filesep 'tack' filesep 'locMax' filesep locMaxFileName{i}];
    load(filename);
    locTM = locMax;
    
    % Crop
    IActin2 = IActin(posInset(i,1):posInset(i,1)+wInset-1,...
        posInset(i,2):posInset(i,2)+wInset-1);
    ITM2 = ITM(posInset(i,1):posInset(i,1)+wInset-1,...
        posInset(i,2):posInset(i,2)+wInset-1);
    C2 = cat(3, IActin2, ITM2, Z2);
    
    locActin1 = locActin(pos(i,1):pos(i,1)+w-1,pos(i,2):pos(i,2)+w-1);
    locActin2 = locActin(posInset(i,1):posInset(i,1)+wInset-1,...
        posInset(i,2):posInset(i,2)+wInset-1);
    
    locTM1 = locTM(pos(i,1):pos(i,1)+w-1,pos(i,2):pos(i,2)+w-1);
    locTM2 = locTM(posInset(i,1):posInset(i,1)+wInset-1,...
        posInset(i,2):posInset(i,2)+wInset-1);
    
    D1 = D(pos(i,1):pos(i,1)+w-1,pos(i,2):pos(i,2)+w-1);
    D2 = D(posInset(i,1):posInset(i,1)+wInset-1,...
        posInset(i,2):posInset(i,2)+wInset-1);
    
    idxActin1 = find(locActin1 ~= 0 & D1 < 5000);
    idxActin2 = find(locActin2 ~= 0 & D2 < 5000);
    idxTM1 = find(locTM1 ~= 0 & D1 < 5000);
    idxTM2 = find(locTM2 ~= 0 & D2 < 5000);
    
    % Display figure
    figure, imshow(C1);
    [y x] = ind2sub(size(locActin1), idxActin1);
    line(x, y,'LineStyle', 'none', 'Marker', '.', 'Color', 'r','MarkerSize',6);
    [y x] = ind2sub(size(locTM1), idxTM1);
    line(x, y,'LineStyle', 'none', 'Marker', '.', 'Color', 'g','MarkerSize',6);
    c = contourc(double(D1), [0; 5000]);
    n = c(2, 1);
    line(c(1, 2:n+2), c(2, 2:n+2), 'Color', 'w', 'Linewidth', .75);
    line(c(1, n+3:end), c(2, n+3:end), 'Color', 'w', 'Linewidth', .75);
    p = posInset(i, :) - pos(i, :);
    line([p(2), p(2) + wInset, p(2) + wInset, p(2), p(2)],...
        [p(1), p(1), p(1) + wInset, p(1) + wInset, p(1)], 'Color', 'w', 'Linewidth', .75);
    print(gcf, '-depsc' , '-painters', [outputDirectory filesep 'fig3_B' ...
        num2str(i) '.eps']);
    
    % Display inset
    figure, imshow(imresize(C2, 8, 'nearest'));
    [y x] = ind2sub(size(locActin2), idxActin2);
    line(8 * x, 8 * y,'LineStyle', 'none', 'Marker', '.', 'Color', 'r','MarkerSize',20);
    [y x] = ind2sub(size(locTM2), idxTM2);
    line(8 * x, 8 * y,'LineStyle', 'none', 'Marker', '.', 'Color', 'g','MarkerSize',20);
    c = contourc(double(D2), [0; 5000]);
    line(c(1, 2:end) * 8, c(2, 2:end) * 8, 'Color', 'w', 'Linewidth', 2);
    
    print(gcf, '-depsc' , '-painters', [outputDirectory filesep 'fig3_B'...
        num2str(i) '_inset.eps']);
end