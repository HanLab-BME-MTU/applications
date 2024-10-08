function mapDescriptives_OneChan(MD, iChan, maxLayer, chanName, chanTitle, figuresDir, varargin)
% mapDescriptives_OneChan Draw descriptive plots of an activity map of
% the specified channel in movieData.
%
% Usage:
%       mapDescriptives_OneChan(MD, 1, 3, 'Actin', 'Actin', ...
%       fullfile(MD.outputDirectory_, 'mapDescriptives'), 'impute', 0, ...
%       'parpoolNum', 4)
%
% Input:
%       MD          - movieData object
%       iChan       - channel index
%                     iChan = 1x indicates the differenced map (X_t =
%                     X_{t-1}) of channel x.
%       maxLayer    - maximum layer to which activity maps are drawn
%       chanName    - a short name for the channel. eg. 'Actin'
%       chanTitle   - a more detailed name for the channel
%                   eg. 'Velocity (nm/sec)'
%       figuresDir  - a directory where plots are saved as png files
%
% Output: png files are saved in the figuresDir
%
% Option:
%       adf         - if true, augmented Dickey-Fuller tests are performed
%                   and ploted. Default is false.
%       figFlag     - if 'on', matlab figures are ploted. Default is 'off'.
%       impute      - if true, moderate missing values are imputed by using
%                   knnimpute.m function. Default is true.
%       parpoolNum  - number of local parallel pool used during permutation. Default is 4.
%       rseed       - input for running rng('default'); rng(rseed). Default
%                   is 'shuffle'. If it is a specific number, the permutation will give
%                   the same result.
%       numPerm     - number of permutation. Default is 1000.
%       WithN       - if true, it uses an alternative windowSampling result
%                   which is obtained by sampleMovieWindowsWithN.m and includes number
%                   of pixels for each windows. Default is false.
%       omittedWindows  
%                   - window index in which activities will be replaced by
%                   NaN. Default is null.
%       subFrames
%                   - specified frames will be only used.      
%       topograph   - if 'off' topographs are not plotted. Default is 'on'.
%
% Updated: 
% J Noh, 2017/11/05, Add topograph for SD and CV with 1%~99% scale bar. Increase FontSize.
% J Noh, 2017/10/11, raw activities can be smoothed. New option is
% 'movingAvgSmoothing'.
% J Noh, 2017/09/25. Now acf is saved.
% J Noh, 2017/08/26. To deal with differenced channels. 
% Jungsik Noh, 2017/05/23
% Jungsik Noh, 2016/10/18


ip = inputParser;
ip.addParameter('adf', 0);
ip.addParameter('figFlag', 'off');
ip.addParameter('impute', true);
ip.addParameter('WithN', false);
ip.addParameter('parpoolNum', 4);
ip.addParameter('rseed', 'shuffle');
ip.addParameter('numPerm', 1000);
ip.addParameter('omittedWindows', []);
ip.addParameter('Folding', false);
ip.addParameter('subFrames', []);
ip.addParameter('topograph', 'on');
ip.addParameter('movingAvgSmoothing', false);


ip.parse(varargin{:});
p = ip.Results;

%figFlag = p.figFlag;


%%  figuresDir setup
% figuresDir = fullfile(outDir, figDirName)           %% input
if ~isdir(figuresDir); mkdir(figuresDir); end

tmptext = ['mapDescriptives_OneChan_', 'inputParser.mat'];
save(fullfile(figuresDir, tmptext), 'p')

%%  getting Maps from channels

disp(chanName)
disp(chanTitle)


[fname0, MDpixelSize_, MDtimeInterval_, wmax, tmax, rawActmap, actmap_outl, imActmap] ...
        = mapOutlierImputation(MD, iChan, maxLayer, 'impute', p.impute, ...
            'WithN', p.WithN, 'omittedWindows', p.omittedWindows, ...
            'Folding', p.Folding, 'subFrames', p.subFrames, 'movingAvgSmoothing', p.movingAvgSmoothing); 
        

        
disp(['== MDpixelSize_: ', num2str(MDpixelSize_), ' =='])
disp(['== MDtimeInterval_: ', num2str(MDtimeInterval_), ' =='])


%%  .txt (export comma delimited files)

for indL = 1:maxLayer
    dlmwrite(fullfile(figuresDir, [fname0, '_', num2str(indL), 'L_actmap_outl.txt']), ...
                    actmap_outl{indL}, 'precision', 8)
    dlmwrite(fullfile(figuresDir, [fname0, '_', num2str(indL), 'L_imActmap.txt']), ...
                    imActmap{indL}, 'precision', 8)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plots 



%%  Raw non-smoothActivityMap prot/act maps

fchanraw = cell(1, maxLayer);
for indL = 1:maxLayer

    inputmap = rawActmap{indL};
    %filteredmap = smoothActivityMap(velmap, 'SmoothParam', smParam, 'UpSample', 1);
    fchanraw{indL} = figure('Visible', p.figFlag);  
    figtmp = imagesc(inputmap);
    title([chanTitle, '-', num2str(indL), 'L'])
    colorbar;colormap(jet) 

    figtmp.AlphaData = 1-isnan(inputmap);
    axis xy;xlabel('Time (s)');ylabel('Window')
    ax = gca;
    curTick = ax.XTick;
    ax.XTickMode = 'manual';
    ax.XTick = curTick+1;
    ax.XTickLabel = (curTick)*MDtimeInterval_;
    
    ax.FontSize = 15;

end


%%
for indL = 1:maxLayer
    saveas(fchanraw{indL}, fullfile(figuresDir, ['/raw', fname0, 'Map_', num2str(indL), 'L.png']), 'png')
end


%%  outl non-smoothActivityMap prot/act maps

fchanraw = cell(1, maxLayer);
for indL = 1:maxLayer

    inputmap = actmap_outl{indL};
    %filteredmap = smoothActivityMap(velmap, 'SmoothParam', smParam, 'UpSample', 1);
    fchanraw{indL} = figure('Visible', p.figFlag);  
    figtmp = imagesc(inputmap);
    title([chanTitle, '-', num2str(indL), 'L'])
    colorbar;colormap(jet) 

    figtmp.AlphaData = 1-isnan(inputmap);
    axis xy;xlabel('Time (s)');ylabel('Window')
    ax = gca;
    curTick = ax.XTick;
    ax.XTickMode = 'manual';
    ax.XTick = curTick+1;
    ax.XTickLabel = (curTick)*MDtimeInterval_;
    
    ax.FontSize = 15;

end


%%
for indL = 1:maxLayer
    saveas(fchanraw{indL}, fullfile(figuresDir, ['outl_', fname0, 'Map_', num2str(indL), 'L.png']), 'png')
end


%%  Run:
%%  BoxPlot per layer

%bpmap = actmap_outl; 
bpmat = nan(wmax*tmax, maxLayer);
for l = 1:maxLayer; 
    bpmat(:, l) = reshape(actmap_outl{l}, [], 1); 
end;
BPlayer = figure('Visible', p.figFlag); 
boxplot(bpmat);
title(chanTitle)
xlabel('Layers')


%%
saveas(BPlayer, fullfile(figuresDir, [fname0, 'BPlayer.png']), 'png')


%%  Boxplots after truncation and outlier detection

velBoxTime = cell(1, maxLayer);
velBoxWin = cell(1, maxLayer);
for indL = 1:maxLayer
    
    inputmap = actmap_outl{indL};
    %
    velBoxTime{indL} = figure('Visible', p.figFlag);
    boxplot(inputmap, 'whisker', Inf);
    title([chanTitle, '-', num2str(indL), 'L'])
    xlabel('Time frame')
    set(gca, 'XTick', 10:10:tmax)
    set(gca, 'XTickLabel', 10:10:tmax)


    velBoxWin{indL} = figure('Visible', p.figFlag); 
    boxplot(inputmap', 'whisker', Inf);
    title([chanTitle, '-', num2str(indL), 'L'])
    xlabel('Window')
    set(gca, 'XTick', 10:10:wmax)
    set(gca, 'XTickLabel', 10:10:wmax)
    
    ax = gca;
    ax.FontSize = 15;

end

%%
for indL = 1:maxLayer
    saveas(velBoxTime{indL}, fullfile(figuresDir, [fname0, 'BoxTime_', num2str(indL), 'L.png']), 'png')
end
pause(0.2)
for indL = 1:maxLayer
    saveas(velBoxWin{indL}, fullfile(figuresDir, [fname0, 'BoxWin_', num2str(indL), 'L.png']), 'png')
end


%%  Histogram 

chanHist = cell(1, maxLayer);
for indL = 1:maxLayer
    
    inputmap = actmap_outl{indL};
    
    sm = summary(inputmap(:));
    title2 = ['m=', sprintf('%0.2f', sm.mean), ' std=', sprintf('%0.2f', sm.std)];
    chanHist{indL} = figure('Visible', p.figFlag);
    histogram(inputmap(:));
    title1 = [chanTitle, '-', num2str(indL), 'L'];
    title({title1, title2})
    
    ax = gca;
    ax.FontSize = 15;

end


%%
for indL = 1:maxLayer
    saveas(chanHist{indL}, fullfile(figuresDir, [fname0, 'Hist_', num2str(indL), 'L.png']), 'png')
end



%% topomap topographMD
if strcmp(p.topograph, 'on')

iWinProc = MD.getProcessIndex('WindowingProcess',1,0);
nBandMax_ = MD.processes_{iWinProc}.nBandMax_;
%
topoMap = nan(wmax, nBandMax_);

for indL = 1:maxLayer
    topoMap(:, indL) = mean(actmap_outl{indL}, 2, 'omitnan');
end

title0 = chanTitle;

topomapFig = topographMD(MD, tmax, 1, topoMap, title0, p.figFlag);
    
    ax = gca;
    ax.FontSize = 15;

%%
saveas(topomapFig, fullfile(figuresDir, ['/topograph_', fname0, '.png']), 'png')

end


%% 2017/11/04 update

histMeans = figure('Visible', p.figFlag);
histogram(topoMap(:), 'BinMethod', 'fd');
title([chanTitle, '-means'])

    ax = gca;
    ax.FontSize = 15;

saveas(histMeans, fullfile(figuresDir, ['/histMeans_', fname0, '.png']), 'png')



%%  smoothActivityMap prot/act maps

smParam = 0.8;

fchan = cell(1, maxLayer);
for indL = 1:maxLayer
    
    inputmap = actmap_outl{indL};

    filteredmap = smoothActivityMap(inputmap, 'SmoothParam', smParam, 'UpSample', 1);
    fchan{indL} = figure('Visible', p.figFlag);
    figtmp = imagesc(filteredmap);
    title([chanTitle, '-', num2str(indL), 'L'])
    colorbar;colormap(jet)

    figtmp.AlphaData = 1-isnan(inputmap);
    axis xy;xlabel('Time (s)');ylabel('Window')
    ax = gca;
    curTick = ax.XTick;
    ax.XTickMode = 'manual';
    ax.XTick = curTick+1;
    ax.XTickLabel = (curTick)*MDtimeInterval_;
    
    ax.FontSize = 15;

end


%% 
for indL = 1:maxLayer
    saveas(fchan{indL}, fullfile(figuresDir, [fname0, 'Map_', num2str(indL), 'L.png']), 'png')
end

 
%%  Means plot

% along time frame

cmeansZ = cell(maxLayer,1);
cmeansZmat = nan(tmax, maxLayer);
for indL = 1:maxLayer
    
    inputmap = actmap_outl{indL};
    
    cmeans = mean(inputmap, 1, 'omitnan');
    cmeansZ{indL} = nanZscore(cmeans);
    cmeansZmat(:, indL) = cmeansZ{indL}';
end
   
    timeFr = 1:tmax;
    timeAxis = (timeFr-1)*MDtimeInterval_;
   
    meansTime = figure('Visible', p.figFlag);
    plot(timeAxis, cmeansZmat);
    title([chanTitle, '-', 'means'])
    xlabel('Time (s)'); ylabel('Z-score')
    legend('1L', '2L', '3L', '4L', '5L', 'Location','northoutside','Orientation','horizontal');
    
    ax = gca;
    ax.FontSize = 15;
    
% along windows

cmeansZ = cell(maxLayer,1);
cmeansZmat = nan(wmax, maxLayer);
for indL = 1:maxLayer
    
    inputmap = actmap_outl{indL};
    
    cmeans = mean(inputmap, 2, 'omitnan');
    cmeansZ{indL} = nanZscore(cmeans);
    cmeansZmat(:, indL) = cmeansZ{indL}';
end
   
    winIndex = 1:wmax;
   
    meansWin = figure('Visible', p.figFlag);
    plot(winIndex, cmeansZmat);
    title([chanTitle, '-', 'means'])
    xlabel('Window'); ylabel('Z-score')
    legend('1L', '2L', '3L', '4L', '5L', 'Location','northoutside','Orientation','horizontal');

    ax = gca;
    ax.FontSize = 15;
    
    
%% 
saveas(meansTime, fullfile(figuresDir, ['/meansTime', fname0, '.png']), 'png')
saveas(meansWin, fullfile(figuresDir, ['/meansWin', fname0, '.png']), 'png')


%%  TS plots for sampled 6 windows

inputmap = actmap_outl{1};
indNotAllNaN = find(~all(isnan(inputmap), 2));
ind0 = round(linspace(1, numel(indNotAllNaN), 6));
winInd = indNotAllNaN(ind0);

tsplots1 = cell(1, maxLayer);
tsplots2 = cell(1, maxLayer);

for indL = 1:maxLayer
    
    inputmap = actmap_outl{indL};

    legend1 = {['win', num2str(winInd(1))], ['win', num2str(winInd(2))], ['win', num2str(winInd(3))]};

    tsplots1{indL} = figure('Visible', p.figFlag);
    plot(timeAxis, inputmap(winInd(1:3), :))
    xlabel('Time (s)'); ylabel(chanName)
    title([chanTitle, ' example1 - ', num2str(indL), 'L'])
    legend(legend1, 'Location','northoutside','Orientation','horizontal')
    
    ax = gca;
    ax.FontSize = 15;

    legend2 = {['win', num2str(winInd(4))], ['win', num2str(winInd(5))], ['win', num2str(winInd(6))]};
    
    tsplots2{indL} = figure('Visible', p.figFlag);
    plot(timeAxis, inputmap(winInd(4:6), :))
    xlabel('Time (s)'); ylabel(chanName)
    title([chanTitle, ' example2 - ', num2str(indL), 'L'])
    legend(legend2, 'Location','northoutside','Orientation','horizontal')
    
    ax = gca;
    ax.FontSize = 15;
end


%% 
for indL = 1:maxLayer
    saveas(tsplots1{indL}, fullfile(figuresDir, ['/tsplots1_', fname0, num2str(indL), 'L.png']), 'png')
    saveas(tsplots2{indL}, fullfile(figuresDir, ['/tsplots2_', fname0, num2str(indL), 'L.png']), 'png')
end



%%  spatial/temporal AutoCorr 1

acmap = cell(1, maxLayer);
corrMat = cell(1, maxLayer);
meanAutocorr = cell(1, maxLayer);
for indL = 1:maxLayer

    mapName = [chanTitle, '-', num2str(indL), 'L'];
    
[acmap{indL}, corrMat{indL}, meanAutocorr{indL}] = ...
    TimeSpaceAutoCorPlot(imActmap{indL}, mapName, MDtimeInterval_, 'figFlag', p.figFlag);
    
end



%%
for indL = 1:maxLayer
    saveas(acmap{indL}, fullfile(figuresDir, ['/acmap', fname0, '_', num2str(indL), 'L.png']), 'png')
    saveas(corrMat{indL}, fullfile(figuresDir, ['/corrMat', fname0, '_', num2str(indL), 'L.png']), 'png')
    saveas(meanAutocorr{indL}, fullfile(figuresDir, ['/meanAutocorr', fname0, '_', num2str(indL), 'L.png']), 'png')
end


%% acCurve = autoCorrCurvePermTest(map, mapName, MDtimeInterval_, numPerm, parpoolNum, rseed)

%numPerm = 1000;
%parpoolNum = 6;
%rseed = 'shuffle';

acCurve = cell(1, maxLayer);
Avg_autocorLayers = cell(1, maxLayer);
for indL = 1:maxLayer
        mapName = [chanTitle, '-', num2str(indL), 'L'];
        
[acCurve{indL}, Avg_autocorLayers{indL}] = autoCorrCurvePermTest(imActmap{indL}, mapName, MDtimeInterval_, ...
                p.numPerm, p.parpoolNum, p.rseed, 'figFlag', p.figFlag);
            
            ax = gca;
            ax.FontSize = 15;
end

%%
for indL=1:maxLayer
saveas(acCurve{indL}, fullfile(figuresDir, ['acCurve_', fname0, '_', num2str(indL), 'L.png']), 'png')
end

%%  09/25/2017
save(fullfile(figuresDir, [fname0, '_Avg_autocorLayers.mat']), 'Avg_autocorLayers')



%%  adftest map

if p.adf == 1

    pvec = cell(1, maxLayer);
    adfMap = cell(1, maxLayer);
    for indL = 1:maxLayer

        mapName = [chanName, '-', num2str(indL), 'L'];
        % 2017.8 p.derivative
        mapForAdf = imActmap{indL}(:, ~all(isnan(imActmap{indL})));
        
        [pvec{indL}, adfMap{indL}] = nanAdfTestMap(mapForAdf, mapName, 0.8);

        ax = gca;
        ax.FontSize = 15;

    end

    % plot
    for indL = 1:maxLayer
        saveas(adfMap{indL}, fullfile(figuresDir, ['/adfMap', fname0, '_', num2str(indL), 'L.png']), 'png')  
    end

end


%%  Coefficient of variation (sd/|mean|)

fvariation = cell(1, maxLayer);
for indL = 1:maxLayer

    map = imActmap{indL};

    % temporal
    stdTemp = std(map, [], 2, 'omitnan');
    meanTemp = mean(map, 2, 'omitnan');
    temporalCV = stdTemp./abs(meanTemp);

    % spatial
    stdSpat = std(map, [], 1, 'omitnan');
    meanSpat = mean(map, 1, 'omitnan');
    spatialCV = stdSpat./abs(meanSpat);

    CVfull = [temporalCV; spatialCV'];
%    stdFull = [stdTemp ; stdSpat'];

    tmax_ = size(map, 2); 
    wmax_ = size(map, 1);

    grChar = cell(wmax_+tmax_, 1);
    for k = 1:wmax_ 
        grChar{k} = 'Temporal';
    end
    for k = 1:tmax_
        grChar{wmax_+k} = 'Spatial';
    end

    fvariation{indL} = figure('Visible', p.figFlag);
    boxplot(CVfull, grChar)
    title([chanName, '-', num2str(indL), 'L'])
    ylabel('Coeff. Variation')
    

end


%%
for indL = 1:maxLayer
    saveas(fvariation{indL}, fullfile(figuresDir, ['/CV', fname0, '_', num2str(indL), 'L.png']), 'png')
end


%%  temporal SD plot

SDlayers = cell(1, maxLayer);
CVlayers = cell(1, maxLayer);

for indL = 1:maxLayer

    map = imActmap{indL};
    % 
    stdTemp = std(map, [], 2, 'omitnan');
    meanTemp = mean(map, 2, 'omitnan');
    temporalCV = stdTemp./abs(meanTemp);

    SDlayers{indL} = stdTemp;
    CVlayers{indL} = temporalCV;
end


%%

fSDlayers = figure('Visible', p.figFlag);

plot(cell2mat(SDlayers))
xlabel('Window'); 
title('Standard deviation')
legend('1L', '2L', '3L', '4L', '5L', 'Location','northoutside','Orientation','horizontal');

ax = gca;
    ax.FontSize = 15;
    

fCVlayers = figure('Visible', p.figFlag);

plot(cell2mat(CVlayers))
xlabel('Window'); 
title('Coeff. Variation')
legend('1L', '2L', '3L', '4L', '5L', 'Location','northoutside','Orientation','horizontal');

ax = gca;
    ax.FontSize = 15;
    

%% save
saveas(fSDlayers, fullfile(figuresDir, ['SDlayers_', fname0, '.png']), 'png')
saveas(fCVlayers, fullfile(figuresDir, ['CVlayers_', fname0, '.png']), 'png')



%% topomap topographMD
if strcmp(p.topograph, 'on')

iWinProc = MD.getProcessIndex('WindowingProcess',1,0);
nBandMax_ = MD.processes_{iWinProc}.nBandMax_;
%
topoMap = nan(wmax, nBandMax_);

for indL = 1:maxLayer
    %topoMap(:, indL) = mean(actmap_outl{indL}, 2, 'omitnan');
    topoMap(:, indL) = SDlayers{indL};
end

title0 = [chanTitle, '-SD'];

topomapFig = topographMD(MD, tmax, 1, topoMap, title0, p.figFlag);

ax = gca;
    ax.FontSize = 15;

%%
saveas(topomapFig, fullfile(figuresDir, ['/topographSD_', fname0, '.png']), 'png')

end



%% 2017/11/04 update

histMeans = figure('Visible', p.figFlag);
histogram(topoMap(:), 'BinMethod', 'fd');
title([chanTitle, '-SDs'])

    ax = gca;
    ax.FontSize = 15;

saveas(histMeans, fullfile(figuresDir, ['/histSDs_', fname0, '.png']), 'png')




%%
if strcmp(p.topograph, 'on')

iWinProc = MD.getProcessIndex('WindowingProcess',1,0);
nBandMax_ = MD.processes_{iWinProc}.nBandMax_;
%
topoMap = nan(wmax, nBandMax_);

for indL = 1:maxLayer
    %topoMap(:, indL) = mean(actmap_outl{indL}, 2, 'omitnan');
    topoMap(:, indL) = CVlayers{indL};
end

title0 = [chanTitle, '-CV'];

topomapFig = topographMD(MD, tmax, 1, topoMap, title0, p.figFlag);

ax = gca;
    ax.FontSize = 15;

%%
saveas(topomapFig, fullfile(figuresDir, ['/topographCV_', fname0, '.png']), 'png')

end


%% 2017/11/04 update

histMeans = figure('Visible', p.figFlag);
histogram(topoMap(:), 'BinMethod', 'fd');
title([chanTitle, '-CVs'])

    ax = gca;
    ax.FontSize = 15;

saveas(histMeans, fullfile(figuresDir, ['/histCVs_', fname0, '.png']), 'png')





%%  checkWindowJump

winTrack = checkWindowJump(MD, 'figFlag', p.figFlag);

%
saveas(winTrack, fullfile(figuresDir, 'window1Trajectory.png'), 'png');


%%
disp('====End of mapDescriptives_OneChan====')
    
end


