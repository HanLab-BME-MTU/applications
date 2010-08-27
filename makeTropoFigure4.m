function makeTropoFigure4(analysisPaths, outputDirectory)

outputDirectory = fullfile(outputDirectory, 'Fig4');

if ~exist(outputDirectory,'dir')
    mkdir(outputDirectory);
end

pathForPanelA = {analysisPaths.ActinTM2{3}, analysisPaths.ActinTM4{1}, analysisPaths.ActinTM5{2}};

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

% Color of inset frame
insetFrameColor = [.798 .057 .259];

for iTM = 1:3
    % Load Movie Data
    fileName = fullfile(pathForPanelA{iTM}, 'movieData.mat');
    if ~exist(fileName, 'file')
        error(['Unable to locate ' fileName]);
    end
    load(fileName);
    
    % Read the list of TMs images (channel 1)
    imagePath = fullfile(movieData.fsmDirectory{1}, 'crop');
    imageFiles = dir([imagePath filesep '*.tif']);
    
    yRange = imagePos(iTM,1):imagePos(iTM,1)+imageSize(1)-1;
    xRange = imagePos(iTM,2):imagePos(iTM,2)+imageSize(2)-1;
    
    %-----------------------------------------------------------------%
    %                                                                 %
    %                          FIGURE 4 PANEL A                       %
    %                                                                 %
    %-----------------------------------------------------------------%
    
    % Read image
    fileName = fullfile(imagePath, imageFiles(iFrames(iTM)).name);
    I = imread(fileName);
    % Read the MPM
    load(fullfile(movieData.fsmDirectory{1}, 'tack', 'mpm.mat'));
    % Crop image
    I = I(yRange,xRange);
    % Crop TM speckles
    idxLoc1 = find(loc1(yRange,xRange) ~= 0);    
    % Crop Actin speckles
    idxLoc2 = find(loc2(yRange,xRange) ~= 0);
    
    % Save
    hFig = figure('Visible', 'off');
    % Set a black frame on the image
    BWima(:,[1,end]) = false;
    BWima([1,end],:) = false;    
    imshow(BWima,[]); hold on;
    [y x] = ind2sub(size(BWima), idxLoc1);    
    line(x, y,'LineStyle', 'none', 'Marker', '.', 'Color', 'g','MarkerSize',12);
    % Draw Actin speckles
    [y x] = ind2sub(size(BWima), idxLoc2);
    line(x, y,'LineStyle', 'none', 'Marker', '.', 'Color', 'r','MarkerSize',12);
    % Draw the inset box
    p = insetPos(iTM, :) - imagePos(iTM, :);
    line([p(2), p(2) + insetSize(2), p(2) + insetSize(2), p(2), p(2)], ...
        [p(1), p(1), p(1) + insetSize(1), p(1) + insetSize(1), p(1)], ...
        'Color', insetFrameColor, 'Linewidth', 2);
    fileName = fullfile(outputDirectory, ['Fig3_B' num2str(iTM) '1.eps']);
    print(hFig, '-depsc' , '-painters', fileName);
    % Close the figure
    close(hFig);
end

% names = fieldnames(analysisPaths);
% paths = cellfun(@(x) analysisPaths.(x), names, 'UniformOutput',false);
% nMovies = cellfun(@numel, paths);
% 
% nBands = 4;
% dLims = [0 .5 1 2 10] * 1000;
% 
% dataB = cell(3,nBands);
% 
% for iTM = 1:numel(paths)
%     % Load Movie Data
%     fileName = [paths{iTM} filesep 'movieData.mat'];
%     if ~exist(fileName, 'file')
%         error(['Unable to locate ' fileName]);
%     end
%     load(fileName);
% 
%     %Verify that the distance transforms have been performed
%     if ~checkMovieBWDist(movieData)
%         error('Distance transforms need to be computed before processing figure 4.');
%     end
%     
%     nrows = movieData.imSize(1);
%     ncols = movieData.imSize(2);
%     nFrames = movieData.labels.nFrames;
%     pixelSize = movieData.pixelSize_nm;
%     timeInterval = movieData.timeInterval_s;
% 
%     % Read the list of distance transforms
%     bwdistPath = movieData.bwdist.directory;
%     bwdistFiles = dir([bwdistPath filesep '*.mat']);
% 
%     % Read distance transforms
%     distToEdge = zeros(nrows,ncols,nFrames);
%     
%     for iFrame = 1:nFrames
%         fileName = fullfile(bwdistPath, bwdistFiles(iFrame).name);
%         tmp = load(fileName);
%         distToEdge(:,:,iFrame) = tmp.distToEdge * pixelSize;
%     end
% 
%     %-----------------------------------------------------------%
%     %                                                           %
%     %                    DATA FOR PANEL B                       %
%     %                                                           %
%     %-----------------------------------------------------------%
% 
%     % Read the MPM
%     load(fullfile(movieData.fsmDirectory{1}, 'tack', 'mpm.mat'));    
% 
%     trackInfos = mpm2trackInfos(MPM,distToEdge,dLims,2);
% 
%     for iBand = 1:nBands
%         dataB{iTM,iBand} = trackInfos{iBand}(:,3) - trackInfos{iBand}(:,2) + 1;
%     end
% end
% 
% %-----------------------------------------------------------------%
% %                                                                 %
% %                          FIGURE 3 PANEL B                       %
% %                                                                 %
% %-----------------------------------------------------------------%
% 
% colors = [
%    0.983333333333333   1.000000000000000   0.800000000000000;
%    0.360000000000000   0.630000000000000   0.900000000000000;
%    0.060000000000000   0.330000000000000   0.600000000000000;
%    0.700000000000000   0.245000000000000   0.245000000000000;
%    0.550000000000000                   0                   0;
%    0.250000000000000                   0                   0; ];
% 
% for iBand = 1:nBands
%     
%     lifeTimeRange = 2:10;
%     
%     hFig = figure('Visible', 'off');
%     set(gca, 'FontName', 'Helvetica', 'FontSize', 20);
%     set(gcf, 'Position', [680 678 560 400], 'PaperPositionMode', 'auto');
%     n = cell2mat(arrayfun(@(iTM) hist(dataB{iTM,iBand}, lifeTimeRange), (1:3)', 'UniformOutput', false));
%     n = bsxfun(@rdivide,n,sum(n,2));
%     h = bar(lifeTimeRange, n','group');
%     
%     set(gca,'XLim',[lifeTimeRange(1)-1, lifeTimeRange(end)+1]);
%     
%     for i=1:numel(h)
%         hC = get(h(i), 'Children');
%         set(hC,'FaceColor', colors(i,:));
%     end
%     
%     set(gca,'XTickLabel',[]);
%     arrayfun(@(x) text(x, -.8/800,[num2str((x-1)*timeInterval) '-' num2str(x*timeInterval) 's'],...
%         'VerticalAlignment', 'Top', 'HorizontalAlignment', 'Right', 'Rotation', 45,...
%         'FontSize', 12), lifeTimeRange(1:end-1));
%     text(lifeTimeRange(end), -.8/800, ['\geq' num2str(lifeTimeRange(end-1) * timeInterval) 's'],...
%         'VerticalAlignment', 'Top', 'HorizontalAlignment', 'Right', 'Rotation', 45,...
%         'FontSize', 12);
%     
%     ylabel('%');
%     
%     text(6, 0.3, [num2str(dLims(iBand)/1000) '-' ...
%         num2str(dLims(iBand+1)/1000) char(181) 'm'], ...
%         'HorizontalAlignment','center', 'FontSize', 20);
%     
%     h = legend({'TM2', 'TM4', 'TM5NM1'});
%     hC = get(h,'Children');
%     hC = hC(1:2:end);
%     for i=1:numel(hC)
%         hCC = get(hC(i),'Children');
%         set(hCC,'FaceColor', colors(numel(hC) - i + 1,:));
%     end
%     legend('boxoff');
%     
%     fileName = [outputDirectory filesep 'Fig4_B' num2str(iBand) '.eps'];
%     print(hFig, '-depsc', fileName);
%     fixEpsFile(fileName);
%     close(hFig);
% end
