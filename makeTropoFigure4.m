function makeTropoFigure4(paths, outputDirectory)

outputDirectory = fullfile(outputDirectory, 'Fig4');

if ~exist(outputDirectory,'dir')
    mkdir(outputDirectory);
end

nBands = 4;
dLims = [0 .5 1 2 10] * 1000;

dataB = cell(3,nBands);

for iTM = 1:numel(paths)
    % Load Movie Data
    fileName = [paths{iTM} filesep 'movieData.mat'];
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

    %-----------------------------------------------------------%
    %                                                           %
    %                    DATA FOR PANEL B                       %
    %                                                           %
    %-----------------------------------------------------------%

    % Read the MPM
    load(fullfile(movieData.fsmDirectory{1}, 'tack', 'mpm.mat'));    

    trackInfos = mpm2trackInfos(MPM,distToEdge,dLims,2);

    for iBand = 1:nBands
        dataB{iTM,iBand} = trackInfos{iBand}(:,3) - trackInfos{iBand}(:,2) + 1;
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

for iBand = 1:nBands
    
    lifeTimeRange = 2:10;
    
    hFig = figure('Visible', 'off');
    set(gca, 'FontName', 'Helvetica', 'FontSize', 20);
    set(gcf, 'Position', [680 678 560 400], 'PaperPositionMode', 'auto');
    n = cell2mat(arrayfun(@(iTM) hist(dataB{iTM,iBand}, lifeTimeRange), (1:3)', 'UniformOutput', false));
    n = bsxfun(@rdivide,n,sum(n,2));
    h = bar(lifeTimeRange, n','group');
    
    set(gca,'XLim',[lifeTimeRange(1)-1, lifeTimeRange(end)+1]);
    
    for i=1:numel(h)
        hC = get(h(i), 'Children');
        set(hC,'FaceColor', colors(i,:));
    end
    
    set(gca,'XTickLabel',[]);
    arrayfun(@(x) text(x, -.8/800,[num2str((x-1)*timeInterval) '-' num2str(x*timeInterval) 's'],...
        'VerticalAlignment', 'Top', 'HorizontalAlignment', 'Right', 'Rotation', 45,...
        'FontSize', 12), lifeTimeRange(1:end-1));
    text(lifeTimeRange(end), -.8/800, ['\geq' num2str(lifeTimeRange(end-1) * timeInterval) 's'],...
        'VerticalAlignment', 'Top', 'HorizontalAlignment', 'Right', 'Rotation', 45,...
        'FontSize', 12);
    
    ylabel('%');
    
    text(6, 0.3, [num2str(dLims(iBand)/1000) '-' ...
        num2str(dLims(iBand+1)/1000) char(181) 'm'], ...
        'HorizontalAlignment','center', 'FontSize', 20);
    
    h = legend({'TM2', 'TM4', 'TM5NM1'});
    hC = get(h,'Children');
    hC = hC(1:2:end);
    for i=1:numel(hC)
        hCC = get(hC(i),'Children');
        set(hCC,'FaceColor', colors(numel(hC) - i + 1,:));
    end
    legend('boxoff');
    
    fileName = [outputDirectory filesep 'Fig4_B' num2str(iBand) '.eps'];
    print(hFig, '-depsc', fileName);
    fixEpsFile(fileName);
    close(hFig);
end
