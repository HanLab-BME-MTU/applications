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
% Location of the crop in Panel A
imagePos = [157 144; 180 67; 30 30];

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
    
    % Extract vectors
    subM = M(:,:,iFrames(iTM) + (-2:2)); %#ok<NODEF>
    
    % Central frame
    tmp = subM(:,:,2);
    
    vectors = tmp(tmp(:,1) >= yRange(1) & tmp(:,1) <= yRange(end) & ...
        tmp(:,2) >= xRange(1) & tmp(:,2) <= xRange(end) & tmp(:,3) ~= 0,:);
    
    xy = vectors(:,1:2);    
    
    for i = [1, 3]
        tmp = subM(:,:,i);
        
        vectors = cat(1, vectors, tmp(tmp(:,1) >= yRange(1) & ...
            tmp(:,1) <= yRange(end) & tmp(:,2) >= xRange(1) & ...
            tmp(:,2) <= xRange(end) & tmp(:,3) ~= 0,:));
    end
    
    Md = vectorFieldAdaptInterp(vectors, xy, 33, [], 'strain');

    % Shift vectors
    Md(:,[1 3]) = Md(:,[1 3]) - yRange(1) + 1;
    Md(:,[2 4]) = Md(:,[2 4]) - xRange(1) + 1;   
    
    % Save
    hFig = figure('Visible', 'off');

    imshow(I,[]); hold on;
    
    Md(:,3:4) = Md(:,1:2) + 10 * (Md(:,3:4) - Md(:,1:2));
    
    quiver(Md(:,2), Md(:,1), Md(:,4) - Md(:,2), Md(:,3) - Md(:,1), 0, 'Color', 'y');
    
%     % Draw the inset box
%     p = insetPos(iTM, :) - imagePos(iTM, :);
%     line([p(2), p(2) + insetSize(2), p(2) + insetSize(2), p(2), p(2)], ...
%         [p(1), p(1), p(1) + insetSize(1), p(1) + insetSize(1), p(1)], ...
%         'Color', insetFrameColor, 'Linewidth', 2);
    fileName = fullfile(outputDirectory, ['Fig4_A' num2str(iTM) '1.eps']);
    print(hFig, '-depsc' , '-painters', fileName);
    % Close the figure
    close(hFig);
end

names = fieldnames(analysisPaths);
paths = cellfun(@(x) analysisPaths.(x), names, 'UniformOutput',false);
nMovies = cellfun(@numel, paths);
 
dLims = [0 .5; .5 1; 1 2; 2 10] * 1000;
nBands = size(dLims,1);

dataB = cell(2,1);
dataB{1} = arrayfun(@(x) cell(x,nBands), nMovies, 'UniformOutput',false);
dataB{2} = arrayfun(@(x) cell(x,nBands), nMovies, 'UniformOutput',false);

for iTM = 1:numel(names)
    
    for iMovie = 1:nMovies(iTM)
        % Load Movie Data% 
        fileName = [paths{iTM}{iMovie} filesep 'movieData.mat'];
        if ~exist(fileName, 'file')
            error(['Unable to locate ' fileName]);
        end
        load(fileName);
        
        %Verify that the distance transforms have been performed
        if ~checkMovieBWDist(movieData)
            error('Distance transforms need to be computed before processing figure 4.');
        end
        
        nrows = movieData.imSize(1);
        ncols = movieData.imSize(2);
        nFrames = movieData.labels.nFrames;
        pixelSize = movieData.pixelSize_nm;
        timeInterval = movieData.timeInterval_s;
        
        % Read the list of distance transforms
        bwdistPath = movieData.bwdist.directory;
        bwdistFiles = dir([bwdistPath filesep '*.mat']);
        
        % Read distance transforms
        distToEdge = zeros(nrows,ncols,nFrames);
        
        for iFrame = 1:nFrames
            fileName = fullfile(bwdistPath, bwdistFiles(iFrame).name);
            tmp = load(fileName);
            distToEdge(:,:,iFrame) = tmp.distToEdge * pixelSize;
        end
        
        % Read the list of label files
        labelPath = movieData.labels.directory;
        labelFiles = dir([labelPath filesep '*.tif']);
        
        % Read labels
        labels = zeros(nrows,ncols,nFrames,'uint16');
        
        for iFrame = 1:nFrames
            fileName = fullfile(labelPath, labelFiles(iFrame).name);
            labels(:,:,iFrame) = imread(fileName);
        end
        
        % Read the protrusion sample
        fileName = fullfile(movieData.protrusion.directory,...
            movieData.protrusion.samples.fileName);
        
        load(fileName);
        
        %---------------------------------------------------------%
        %                                                         %
        %                      DATA FOR PANEL B                   %
        %                                                         %
        %---------------------------------------------------------%
        
        for iProtState = 1:2
            
            % Read the MPM
            load(fullfile(movieData.fsmDirectory{1}, 'tack', 'mpm.mat'));
            
            % set distToEdge to +inf when cell is not in iProtState
            distToEdge2 = zeros(size(distToEdge));
            
            for iFrame = 1:nFrames-1
                idx = find(protrusionSamples.states(:,iFrame) ~= iProtState+1);
                tmp = distToEdge(:,:,iFrame);
                tmp(ismember(labels(:,:,iFrame), idx)) = +inf;
                distToEdge2(:,:,iFrame) = tmp;
            end
            
            distToEdge2(:,:,end) = +inf;
            
            trackInfos = mpm2trackInfos(MPM,distToEdge2,dLims,2);
            
            for iBand = 1:nBands
                dataB{iProtState,iTM,iMovie,iBand} = trackInfos{iBand}(:,3) - trackInfos{iBand}(:,2) + 1;
            end
        end
    end
end

colors = [
    0.290000000000000   0.617000000000000   0.547000000000000;
    0.360000000000000   0.630000000000000   0.900000000000000;
    0.060000000000000   0.330000000000000   0.600000000000000;
    0.000000000000000   0.000000000000000   0.000000000000000
    0.000000000000000   0.000000000000000   0.000000000000000
    0.000000000000000   0.000000000000000   0.000000000000000];

%-----------------------------------------------------------------%
%                                                                 %
%                          FIGURE 3 PANEL B                       %
%                                                                 %
%-----------------------------------------------------------------%

for iProtState = 1:2
    for iBand = 1:nBands
        
        lifeTimeRange = 2:10;
        
        hFig = figure('Visible', 'off');
        set(gca, 'FontName', 'Helvetica', 'FontSize', 20);
        set(gcf, 'Position', [680 678 560 400], 'PaperPositionMode', 'auto');
        
        n = arrayfun(@(x) cell(x,nBands), nMovies, 'UniformOutput',false);
        
        for iTM = 1:numel(names)
            for iMovie = 1:nMovies(iTM)
                n{iTM,iMovie} = hist(dataB{iProtState,iTM,iMovie,iBand}, lifeTimeRange);
            end
        end
        
        mu = zeros(numel(names), numel(lifeTimeRange));
        sigma = zeros(size(mu));
        
        for iTM = 1:numel(names)
            tmp = vertcat(n{iTM,:});
            tmp = tmp ./ sum(tmp(:));
            mu(iTM,:) = sum(tmp);
            sigma(iTM,:) = std(tmp);
        end
        
        mu = cumsum(mu,2);
        
        hold on;
        for iTM = 1:numel(names)
            line(lifeTimeRange, mu(iTM,:), 'Color', colors(iTM,:), 'LineWidth', 1.5);
        end
        for iTM = 1:numel(names)
            errorbar(lifeTimeRange,mu(iTM,:),sigma(iTM,:),'.', 'Color', colors(iTM,:));
        end
        
        set(gca,'XLim',[lifeTimeRange(1)-1, lifeTimeRange(end)+1]);
        set(gca,'YLim',[-1/80 1]);
        
        set(gca,'XTickLabel',[]);
        arrayfun(@(x) text(x, -.8/800,[num2str((x-1)*timeInterval) '-' num2str(x*timeInterval) 's'],...
            'VerticalAlignment', 'Top', 'HorizontalAlignment', 'Right', 'Rotation', 45,...
            'FontSize', 12,'Parent',gca), lifeTimeRange(1:end-1));
        text(lifeTimeRange(end), -.8/800, ['\geq' num2str(lifeTimeRange(end-1) * timeInterval) 's'],...
            'VerticalAlignment', 'Top', 'HorizontalAlignment', 'Right', 'Rotation', 45,...
            'FontSize', 12,'Parent',gca);
        
        hold off;
        
        legend({'TM2', 'TM4', 'TM5NM1'},'Location', 'Best');
        legend('boxoff');
        
        fileName = [outputDirectory filesep 'Fig4_B' num2str(iBand) num2str(iProtState) '.eps'];
        print(hFig, '-depsc', fileName);
        fixEpsFile(fileName);
        close(hFig);
    end
end
