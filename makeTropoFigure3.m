function makeTropoFigure3(paths, outputDirectory)

outputDirectory = fullfile(outputDirectory, 'Fig3');

if ~exist(outputDirectory,'dir')
    mkdir(outputDirectory);
end

nMovies = numel(paths);

% Chosen frame to display in Panel A
iFrames = [16, 30, 42];
% Size of images in Panel A
imageSize = [390 390];
% Size of the inset in Panel A
insetSize = [180 180];
% Location of the crop in Panel A
imagePos = [157 144; 180 67; 30 30];
% Location of the inset in Panel A
insetPos = [304,267; 275,167; 233 134];

dataC = cell(nMovies,2);
dataD1 = cell(nMovies,1);
dataD2 = cell(nMovies,1);

maxDist = 5000;

for iTM = 1:numel(paths)
    % Load Movie Data
    fileName = [paths{iTM} filesep 'movieData.mat'];
    if ~exist(fileName, 'file')
        error(['Unable to locate ' fileName]);
    end
    load(fileName);

    %Verify that the labeling has been performed
    if ~checkMovieLabels(movieData)
        error('Must label movie before computing figure 3.');
    end

    nrows = movieData.imSize(1);
    ncols = movieData.imSize(2);
    nFrames = movieData.labels.nFrames;
    pixelSize = movieData.pixelSize_nm;
    
    % Read the list of mask files
    maskPath = movieData.masks.directory;
    maskFiles = dir([maskPath filesep '*.tif']);
    
    % Read the list of label files
    labelPath = movieData.labels.directory;
    labelFiles = dir([labelPath filesep '*.tif']);

    % Read the list of TMs speckles (channel 1)
    s1Path = fullfile(movieData.fsmDirectory{1}, 'tack', 'locMax');
    s1Files = dir([s1Path filesep '*.mat']);

    % Read the list of Actin speckles (channel 2)
    s2Path = fullfile(movieData.fsmDirectory{2}, 'tack', 'locMax');
    s2Files = dir([s2Path filesep '*.mat']);

    % Read the list of TMs images (channel 1)
    image1Path = fullfile(movieData.fsmDirectory{1}, 'crop');
    image1Files = dir([image1Path filesep '*.tif']);
    
    % Read the list of Actin images (channel 2)
    image2Path = fullfile(movieData.fsmDirectory{2}, 'crop');
    image2Files = dir([image2Path filesep '*.tif']);
    
    % Read the list of distance transforms
    bwdistPath = movieData.bwdist.directory;
    bwdistFiles = dir([bwdistPath filesep '*.mat']);

    % Load activity map
    fileName = fullfile(movieData.protrusion.directory,...
        movieData.protrusion.samples.fileName);
    if ~exist(fileName, 'file')
        error(['Unable to locate ' fileName]);
    end
    load(fileName);
    v = sort(protrusionSamples.averageMagnitude(:));
    val = v(ceil(.01 * numel(v)));
    protMask = protrusionSamples.averageNormalComponent > val;
    retMask = protrusionSamples.averageNormalComponent < -val;
    
    yRange = imagePos(iTM,1):imagePos(iTM,1)+imageSize(1)-1;
    xRange = imagePos(iTM,2):imagePos(iTM,2)+imageSize(2)-1;
    
    yRangeInset = insetPos(iTM,1):insetPos(iTM,1)+insetSize(1)-1;
    xRangeInset = insetPos(iTM,2):insetPos(iTM,2)+insetSize(2)-1;
    
    %-----------------------------------------------------------------%
    %                                                                 %
    %                          FIGURE 3 PANEL A                       %
    %                                                                 %
    %-----------------------------------------------------------------%
    
    % TM (Panel A, row 1)
    
    % Read image
    fileName = fullfile(image1Path, image1Files(iFrames(iTM)).name);
    I = imread(fileName);
    % Crop image
    I1 = double(I(yRange,xRange));
    I1 = (I1 - min(I1(:))) / max(I1(:));
    % Save
    fileName = fullfile(outputDirectory, ['Fig3_A' num2str(iTM) '1.tif']);
    imwrite(I1,fileName);
    
    % Actin (Panel A, row 2)
    
    % Read image
    fileName = fullfile(image2Path, image2Files(iFrames(iTM)).name);
    I = imread(fileName);
    % Crop image
    I2 = double(I(yRange,xRange));
    I2 = (I2 - min(I2(:))) / max(I2(:));
    % Save
    fileName = fullfile(outputDirectory, ['Fig3_A' num2str(iTM) '2.tif']);
    imwrite(uint8(I2 * 255),fileName);

    % Merge (Panel A, row 3)
    C = cat(3,I2, I1, zeros(size(I1)));
    % Save
    fileName = fullfile(outputDirectory, ['Fig3_A' num2str(iTM) '3.tif']);
    imwrite(C,fileName);
    
    %-----------------------------------------------------------------%
    %                                                                 %
    %                          FIGURE 3 PANEL B                       %
    %                                                                 %
    %-----------------------------------------------------------------%

    % Mask + TM speckle + Actin speckle + inset frame (Panel B, row 1)
    
    % Read mask
    fileName = fullfile(maskPath, maskFiles(iFrames(iTM)).name);
    BW = imread(fileName);
    % Load TM speckles
    fileName = fullfile(s1Path, s1Files(iFrames(iTM)).name);
    load(fileName);
    loc1 = locMax;
    % Load Actin speckles
    fileName = fullfile(s2Path, s2Files(iFrames(iTM)).name);
    load(fileName);
    loc2 = locMax;
    % Crop mask
    BWima = BW(yRange,xRange);
    % Crop TM speckles
    idxLoc1 = find(loc1(yRange,xRange) ~= 0);    
    % Crop Actin speckles
    idxLoc2 = find(loc2(yRange,xRange) ~= 0);
    
    % Save
    hFig = figure('Visible', 'off');
    imshow(BWima,[]); hold on;
    [y x] = ind2sub(size(BWima), idxLoc1);    
    line(x, y,'LineStyle', 'none', 'Marker', '.', 'Color', 'g','MarkerSize',6);
    % Draw Actin speckles
    [y x] = ind2sub(size(BWima), idxLoc2);
    line(x, y,'LineStyle', 'none', 'Marker', '.', 'Color', 'r','MarkerSize',6);
    % Draw the inset box
    p = insetPos(iTM, :) - imagePos(iTM, :);
    line([p(2), p(2) + insetSize(2), p(2) + insetSize(2), p(2), p(2)], ...
        [p(1), p(1), p(1) + insetSize(1), p(1) + insetSize(1), p(1)], ...
        'Color', 'y', 'Linewidth', 2);
    fileName = fullfile(outputDirectory, ['Fig3_B' num2str(iTM) '1.eps']);
    print(hFig, '-depsc' , '-painters', fileName);
    % Close the figure
    close(hFig);
       
    % Inset of Fig3 B.3 (Panel B, row 2)
    
    % Crop image
    Cinset = C(yRangeInset - yRange(1) + 1, xRangeInset - xRange(1) + 1,:);
    fileName = fullfile(outputDirectory, ['Fig3_B' num2str(iTM) '2.tif']);
    imwrite(Cinset,fileName);
    
    % Inset of Fig3 B.4 (Panel B, row 3)

    % Crop mask
    BWinset = BW(yRangeInset, xRangeInset);
    % Crop TM speckles
    idxLoc1 = find(loc1(yRangeInset,xRangeInset) ~= 0);    
    % Crop Actin speckles
    idxLoc2 = find(loc2(yRangeInset,xRangeInset) ~= 0);
    
    % Save
    hFig = figure('Visible', 'off');
    % Set a black frame on the image
    BWinset(:,[1,end]) = false;
    BWinset([1,end],:) = false;
    imshow(BWinset,[]); hold on;
    [y x] = ind2sub(size(BWinset), idxLoc1);    
    line(x, y,'LineStyle', 'none', 'Marker', '.', 'Color', 'g','MarkerSize',6);
    % Drw Actin speckles
    [y x] = ind2sub(size(BWinset), idxLoc2);
    line(x, y,'LineStyle', 'none', 'Marker', '.', 'Color', 'r','MarkerSize',6);
    fileName = fullfile(outputDirectory, ['Fig3_B' num2str(iTM) '3.eps']);
    print(hFig, '-depsc' , '-painters', fileName);
    % Close the figure
    close(hFig);
    
    %-----------------------------------------------------------%
    %                                                           %
    %                    DATA FOR PANEL C                       %
    %                                                           %
    %-----------------------------------------------------------%

    % Read distance transforms
    distToEdge = zeros(nrows,ncols,nFrames);
    
    for iFrame = 1:nFrames
        fileName = fullfile(bwdistPath, bwdistFiles(iFrame).name);
        tmp = load(fileName);
        distToEdge(:,:,iFrame) = tmp.distToEdge * pixelSize;
    end
    
    for iFrame = 1:nFrames-1    
        % Load label
        L = imread(fullfile(labelPath, labelFiles(iFrame).name));
        
        % Load TM speckles (channel 1)
        load(fullfile(s1Path, s1Files(iFrame).name));
        locMax1 = locMax;
        
        % Load Actin speckles (channel 2)
        load(fullfile(s2Path, s2Files(iFrame).name));
        locMax2 = locMax;
        
        % Convert the distance transform in nanometers
        distToEdgeFrame = distToEdge(:,:,iFrame);
    
        % Restrinct locMax to the first maxDist nanometer
        outsideBand = distToEdgeFrame > maxDist;
        locMax1(outsideBand) = 0;
        locMax2(outsideBand) = 0;
                
        Lprot = ismember(L, find(protMask(:, iFrame) == 1));
        
        idxS1 = locMax1 ~= 0 & Lprot;
        idxS2 = locMax2 ~= 0 & Lprot;

        if any(idxS1(:)) && any(idxS2(:))
            dataC{iTM,1} = cat(1, dataC{iTM,1}, mean(distToEdgeFrame(idxS1)) - ...
                mean(distToEdgeFrame(idxS2)));
        end
        
        Lret = ismember(L, find(retMask(:, iFrame) == 1));
        
        idxS1 = locMax1 ~= 0 & Lret;
        idxS2 = locMax2 ~= 0 & Lret;
        
        if any(idxS1(:)) && any(idxS2(:))
            dataC{iTM,2} = cat(1, dataC{iTM,2}, mean(distToEdgeFrame(idxS1)) - ...
                mean(distToEdgeFrame(idxS2)));
        end
    end
    
    %-----------------------------------------------------------%
    %                                                           %
    %                    DATA FOR PANEL D                       %
    %                                                           %
    %-----------------------------------------------------------%

    % Load density scores
    load(fullfile(movieData.density.directory, movieData.density.channelDirectory{1}, 'densityScores.mat'));
    
    distFromEdge = vertcat(densityScores(:).distFromEdge) * pixelSize; %#ok<NODEF>
    averageDensity = vertcat(densityScores(:).averageDensity) * (pixelSize/1000)^-2; % in um-2
    protrusionState = vertcat(densityScores(:).protrusionState);
    %protrusionPersistence = vertcat(densityScores(:).protrusionPersistence);
    %protrusionSpeed = abs(vertcat(densityScores(:).protrusionSpeed) * pixelSize / timeInterval);

    %
    % Data for Panel D1
    %
    
    maxDistFromEdge = min(15000,max(distFromEdge));    
    dist = 0:500:maxDistFromEdge;
    
    dataD1{iTM} = zeros(numel(dist)-1,1);

    for i = 1:numel(dist)-1
        dataD1{iTM}(i) = ...
            mean(averageDensity(distFromEdge > dist(i) & distFromEdge <= dist(i+1) & protrusionState == 2));
    end
    
    %
    % Data for Panel D2
    %

    dataD2{iTM} = zeros(numel(dist)-1,1);
    
    for i = 1:numel(dist)-1
        dataD2{iTM}(i) = ...
            mean(averageDensity(distFromEdge > dist(i) & distFromEdge <= dist(i+1) & protrusionState == 3));
    end
end

%-----------------------------------------------------------------%
%                                                                 %
%                          FIGURE 3 PANEL C                       %
%                                                                 %
%-----------------------------------------------------------------% 

colors = [
   0.983333333333333   1.000000000000000   0.800000000000000;
   0.360000000000000   0.630000000000000   0.900000000000000;
   0.060000000000000   0.330000000000000   0.600000000000000;
   0.700000000000000   0.245000000000000   0.245000000000000;
   0.550000000000000                   0                   0;
   0.250000000000000                   0                   0; ];


hFig = figure('Visible', 'off');
set(gca, 'FontName', 'Helvetica', 'FontSize', 18);
set(gcf, 'Position', [680 678 560 400], 'PaperPositionMode', 'auto');
mu = cellfun(@mean, dataC);
sigma = cellfun(@std, dataC);
h = bar(gca, mu', 'group'); hold on;

XTicks = zeros(size(mu));

for i=1:numel(h)
    hC = get(h(i), 'Children');
    set(hC,'FaceColor', colors(i,:));
    XData = get(hC, 'XData');
    XTicks(i,1:2) = .5 * (XData(1,1:2) + XData(3,1:2));
end

errorbar(XTicks(:),mu(:),sigma(:),'xk'); hold off;

set(gca,'XTickLabel',{'Protrusion', 'Retraction'})
ylabel('Distance to Actin Front (nm)');

h = legend({'TM2', 'TM4', 'TM5NM1'}); legend('boxoff');
hC = get(h,'Children');
hC = hC(1:2:end);
for i=1:numel(hC)
    hCC = get(hC(i),'Children');
    set(hCC,'FaceColor', colors(numel(hC) - i + 1,:));
end

fileName = [outputDirectory filesep 'Fig3_C.eps'];
print(hFig, '-depsc', fileName);
fixEpsFile(fileName);
close(hFig);
   
%-----------------------------------------------------------------%
%                                                                 %
%                          FIGURE 3 PANEL D                       %
%                                                                 %
%-----------------------------------------------------------------% 

%
% Panel D1
%

% convert dist in microns:
dist = dist / 1000;

hFig = figure('Visible', 'off');
set(gca,'XTick',dist(1:4:end-1));
set(gca, 'FontName', 'Helvetica', 'FontSize', 18);
set(gcf, 'Position', [680 678 650 450], 'PaperPositionMode', 'auto');

averageDensityTotal = horzcat(dataD1{:});

% plot x axis in um
h = line(dist(1:end-1), averageDensityTotal,'LineWidth',2);

for i=1:numel(h)
    set(h(i),'Color', colors(i,:));
end

h = legend({'TM2', 'TM4', 'TM5NM1'}); legend('boxoff');
hC = get(h,'Children');
hC = hC(1:2:end);

for i=1:numel(hC)
    hCC = get(hC(i),'Children');
    set(hCC,'FaceColor', colors(numel(hC) - i + 1,:));
end

xlabel(['Distance away from cell edge during protrusion (' char(181) 'm)']);
ylabel(['Speckle Density (' char(181) 'm^{-2})']);

fileName = [outputDirectory filesep 'Fig3_D1.eps'];
print(hFig, '-depsc', fileName);
fixEpsFile(fileName);
close(hFig);

%
% Panel D2
%

hFig = figure('Visible', 'off');
set(gca,'XTick',dist(1:4:end-1));
set(gca, 'FontName', 'Helvetica', 'FontSize', 18);
set(gcf, 'Position', [680 678 650 450], 'PaperPositionMode', 'auto');

averageDensityTotal = horzcat(dataD2{:});

% plot x axis in um
h = line(dist(1:end-1), averageDensityTotal,'LineWidth',2);

for i=1:numel(h)
    set(h(i),'Color', colors(i,:));
end

h = legend({'TM2', 'TM4', 'TM5NM1'}); legend('boxoff');
hC = get(h,'Children');
hC = hC(1:2:end);

for i=1:numel(hC)
    hCC = get(hC(i),'Children');
    set(hCC,'FaceColor', colors(numel(hC) - i + 1,:));
end

xlabel(['Distance away from cell edge during retraction (' char(181) 'm)']);
ylabel(['Speckle Density (' char(181) 'm^{-2})']);

fileName = [outputDirectory filesep 'Fig3_D2.eps'];
print(hFig, '-depsc', fileName);
fixEpsFile(fileName);
close(hFig);
