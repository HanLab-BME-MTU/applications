function makeTropoFigure4(paths, outputDirectory)

% Chosen frame to display in Panel A
iFrames = [16, 74, 42];
% Size of images in Panel A
imageSize = [463 496; 567 624; 419 502];
% Size of the inset in Panel A
insetSize = [180 180; 180 180; 180 180];
% Location of the crop in Panel A
imagePos = [106 142; 1 1; 30 30];
% Location of the inset in Panel A
insetPos = [304,267; 154,161; 233 134];

dataB = cell(3,2);
dataC = cell(3,1);

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

    nFrames = movieData.labels.nFrames;
    pixelSize = movieData.pixelSize_nm;
    timeInterval = movieData.timeInterval_s;
    
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
    
    yRange = imagePos(iTM,1):imagePos(iTM,1)+imageSize(iTM,1)-1;
    xRange = imagePos(iTM,2):imagePos(iTM,2)+imageSize(iTM,2)-1;
    
    yRangeInset = insetPos(iTM,1):insetPos(iTM,1)+insetSize(iTM,1)-1;
    xRangeInset = insetPos(iTM,2):insetPos(iTM,2)+insetSize(iTM,2)-1;
    
    %-----------------------------------------------------------------%
    %                                                                 %
    %                          FIGURE 3 PANEL A                       %
    %                                                                 %
    %-----------------------------------------------------------------%
    
    % TM (Panel A, column 1)
    
    % Read image
    fileName = fullfile(image1Path, image1Files(iFrames(iTM)).name);
    I = imread(fileName);
    % Crop image
    I1 = I(yRange,xRange);
    % Save
    hFig = figure('Visible', 'off');
    imshow(I1, []);
    fileName = fullfile(outputDirectory, ['Fig3_A' num2str(iTM) '1.eps']);
    print(hFig, '-depsc' , '-painters', fileName);
    % Close the figure
    close(hFig);
    
    % Actin (Panel A, column 2)
    
    % Read image
    fileName = fullfile(image2Path, image2Files(iFrames(iTM)).name);
    I = imread(fileName);
    % Crop image
    I2 = I(yRange,xRange);
    % Save
    hFig = figure('Visible', 'off');
    imshow(I2, []);
    fileName = fullfile(outputDirectory, ['Fig3_A' num2str(iTM) '2.eps']);
    print(hFig, '-depsc' , '-painters', fileName);
    % Close the figure
    close(hFig);

    % Mask + TM speckle + Actin speckle + inset frame (Panel A, column 3)
    
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
    line([p(2), p(2) + insetSize(iTM,2), p(2) + insetSize(iTM,2), p(2), p(2)], ...
        [p(1), p(1), p(1) + insetSize(iTM,1), p(1) + insetSize(iTM,1), p(1)], ...
        'Color', [.3 .3 .3], 'Linewidth', 1);
    fileName = fullfile(outputDirectory, ['Fig3_A' num2str(iTM) '3.eps']);
    print(hFig, '-depsc' , '-painters', fileName);
    % Close the figure
    close(hFig);
       
    % Inset of Fig4 A.3 (Panel A, column 4)

    % Crop mask
    BWinset = BW(yRangeInset, xRangeInset);
    % Crop TM speckles
    idxLoc1 = find(loc1(yRangeInset,xRangeInset) ~= 0);    
    % Crop Actin speckles
    idxLoc2 = find(loc2(yRangeInset,xRangeInset) ~= 0);
    
    % Save
    hFig = figure('Visible', 'off');
    imshow(BWinset,[]); hold on;
    [y x] = ind2sub(size(BWinset), idxLoc1);    
    line(x, y,'LineStyle', 'none', 'Marker', '.', 'Color', 'g','MarkerSize',6);
    % Drw Actin speckles
    [y x] = ind2sub(size(BWinset), idxLoc2);
    line(x, y,'LineStyle', 'none', 'Marker', '.', 'Color', 'r','MarkerSize',6);
    fileName = fullfile(outputDirectory, ['Fig3_A' num2str(iTM) '4.eps']);
    print(hFig, '-depsc' , '-painters', fileName);
    % Close the figure
    close(hFig);
    
    %-----------------------------------------------------------%
    %                                                           %
    %                    DATA FOR PANEL B & C                   %
    %                                                           %
    %-----------------------------------------------------------%

    % Read the MPM
    load(fullfile(movieData.fsmDirectory{1}, 'tack', 'mpm.mat'));    
    nrows = size(MPM,1); %#ok<NODEF>   
    trackMask = MPM(:,1:2:end) ~= 0;
    accu = zeros(nrows,1);
    inBand = false(nrows,1);    
    trackID = (1:nrows)';
    
    
    for iFrame = 1:nFrames-1
        %% Data For Panel B
        
        % Load label
        L = imread(fullfile(labelPath, labelFiles(iFrame).name));
        
        % Load TM speckles (channel 1)
        load(fullfile(s1Path, s1Files(iFrame).name));
        locMax1 = locMax;
        
        % Load Actin speckles (channel 2)
        load(fullfile(s2Path, s2Files(iFrame).name));
        locMax2 = locMax;
        
        % Read the distance transform
        fileName = fullfile(bwdistPath, bwdistFiles(iFrame).name);
        load(fileName);
        distToEdge = distToEdge * pixelSize;
    
        % Restrinct locMax to the first maxDist nanometer
        validDist = distToEdge > maxDist;
        locMax1(validDist) = 0;
        locMax2(validDist) = 0;
                
        Lprot = ismember(L, find(protMask(:, iFrame) == 1));
        
        idxS1 = locMax1 ~= 0 & Lprot;
        idxS2 = locMax2 ~= 0 & Lprot;

        if any(idxS1(:)) && any(idxS2(:))
            dataB{iTM,1} = cat(1, dataB{iTM,1}, mean(distToEdge(idxS1)) - mean(distToEdge(idxS2)));
        end
        
        Lret = ismember(L, find(retMask(:, iFrame) == 1));
        
        idxS1 = locMax1 ~= 0 & Lret;
        idxS2 = locMax2 ~= 0 & Lret;
        
        if any(idxS1(:)) && any(idxS2(:))
            dataB{iTM,2} = cat(1, dataB{iTM,2}, mean(distToEdge(idxS1)) - mean(distToEdge(idxS2)));
        end
        
        %% Data for panel C
        
        idxLive = trackID(trackMask(:,iFrame));
        idxDead = trackID(~trackMask(:,iFrame));
        
        % Accumulate lifetime
        accu(idxLive) = accu(idxLive) + 1;
            
        % Check whether idxLive points are within the maxDist band. A track
        % is considered to be within the band (inBand == true) if it has
        % been within the band at least during 1 frame.
        idxLivePoints = sub2ind(size(L), MPM(idxLive,2*iFrame-1), MPM(idxLive,2*iFrame));
        inBand(idxLive) = inBand(idxLive) | locMax1(idxLivePoints) ~= 0;
        
        % Keep
        dataC{iTM} = vertcat(dataC{iTM}, accu(~trackMask(:,iFrame) & accu > 0 & inBand));
        
        % Reset
        accu(idxDead) = 0;
        inBand(idxDead) = false;
    end    
end
 
%-----------------------------------------------------------------%
%                                                                 %
%                          FIGURE 3 PANEL B                       %
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
mu = cellfun(@mean, dataB);
sigma = cellfun(@std, dataB);
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

fileName = [outputDirectory filesep 'Fig3_B.eps'];
print(hFig, '-depsc', fileName);
fixEpsFile(fileName);
close(hFig);
   
%-----------------------------------------------------------------%
%                                                                 %
%                          FIGURE 3 PANEL C                       %
%                                                                 %
%-----------------------------------------------------------------%

lifeTimeRange = 1:10;

hFig = figure('Visible', 'off');
set(gca, 'FontName', 'Helvetica', 'FontSize', 20);
set(gcf, 'Position', [680 678 560 400], 'PaperPositionMode', 'auto');
n = cell2mat(arrayfun(@(i) hist(dataC{i}, lifeTimeRange), (1:3)', 'UniformOutput', false));
n = bsxfun(@rdivide,n,sum(n,2));
h = bar(lifeTimeRange, n','group');

set(gca,'XLim',[lifeTimeRange(1)-1, lifeTimeRange(end)+1]);

for i=1:numel(h)
    hC = get(h(i), 'Children');
    set(hC,'FaceColor', colors(i,:));
end

xTickLabels = arrayfun(@int2str,lifeTimeRange,'UniformOutput',false);
xTickLabels{end} = [xTickLabels{end} '+'];
set(gca,'XTickLabel',xTickLabels);

xlabel('# frame');
ylabel('%');
h = legend({'TM2', 'TM4', 'TM5NM1'});
hC = get(h,'Children');
hC = hC(1:2:end);
for i=1:numel(hC)
    hCC = get(hC(i),'Children');
    set(hCC,'FaceColor', colors(numel(hC) - i + 1,:));
end
legend('boxoff');

fileName = [outputDirectory filesep 'Fig3_C.eps'];
print(hFig, '-depsc', fileName);
fixEpsFile(fileName);
close(hFig);
