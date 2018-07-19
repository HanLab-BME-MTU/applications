function mapDescriptives_Vel(MD, figuresDir, varargin)
% mapDescriptives_Vel Draw descriptive plots for the protrusion map of MD.
%
% Usage:
%       mapDescriptives_Vel(MD, figuresDir, 'adf', 1, 'impute', 1, 'parpoolNum', 4)
%
% Input:
%       MD          - movieData object
%       figuresDir  - a directory where plots are saved as png files
%
% Output:           png files are saved in the figuresDir
%
% Option:
%       adf         - if true, augmented Dickey-Fuller tests are performed
%                   and ploted.
%       figFlag     - if 'on', matlab figures are ploted. Default is 'off'.
%       impute      - if true, moderate missing values are imputed by using
%                   knnimpute.m function.
%       parpoolNum  - number of local parallel pool used during permutation. Default is 4.
%       rseed       - input for running rng('default'); rng(rseed). Default
%                   is 'shuffle'. If it is a specific number, the permutation will give
%                   the same result.
%       numPerm     - number of permutation. Default is 1000.
%       omittedWindows  
%                   - window index in which activities will be replaced by
%                   NaN. Default is null.
%       subFrames
%                   - specified frames will be only used.  
%       derivative  - if true, it computes differenced velocities per
%                   second.
%       topograph   - if 'off' topographs are not plotted. Default is 'on'.
%
% Updated: J Noh, 2017/10/11, raw activities can be smoothed. New option is
% 'movingAvgSmoothing'.
% J Noh, 2017/08/26. To deal with differenced channels. 
% Jungsik Noh, 2017/05/23
% Jungsik Noh, 2016/10/04


ip = inputParser;
ip.addParameter('adf', 1);
ip.addParameter('figFlag', 'off');
ip.addParameter('impute', true);
ip.addParameter('parpoolNum', 4);
ip.addParameter('rseed', 'shuffle');
ip.addParameter('numPerm', 1000);
ip.addParameter('omittedWindows', []);
ip.addParameter('Folding', false);
ip.addParameter('subFrames', []);
ip.addParameter('derivative', false);
ip.addParameter('topograph', 'on');
ip.addParameter('movingAvgSmoothing', false);

ip.parse(varargin{:});
p = ip.Results;

figFlag = p.figFlag;


%%  figuresDir setup
% figuresDir = fullfile(outDir, figDirName)           %% input
if ~isdir(figuresDir); mkdir(figuresDir); end

tmptext = ['mapDescriptives_Vel_', 'inputParser.mat'];
save(fullfile(figuresDir, tmptext), 'p')

%%  getting Maps from channels

switch p.derivative
    case false
        chan0Title = 'Velocity (nm/sec)';
        chan0Name = 'Vel';
        iChan = 0;
    case true
        chan0Title = 'D.Velocity (nm/sec^2)';
        chan0Name = 'D.Vel';
        iChan = 10;
end

%  fnameChan = double(p.derivative)*10 + iChan;

        disp(chan0Name)
        disp(chan0Title)
        disp(iChan);
        maxLayer = 1;

[~, MDpixelSize_, MDtimeInterval_, wmax, tmax, rawActmap, actmap_outl, imActmap] ...
        = mapOutlierImputation(MD, iChan, maxLayer, 'impute', p.impute, ...
            'omittedWindows', p.omittedWindows, 'Folding', p.Folding, ...
            'subFrames', p.subFrames, 'movingAvgSmoothing', p.movingAvgSmoothing); 
        
% ..st layer

    velmap = rawActmap{1};
    velmap_outl = actmap_outl{1};
    imvelocitymap = imActmap{1}(:, 2:tmax);     %  Imputation (used as an input of computations)
                                                    % Note 2:tmax    

if isempty(MDpixelSize_); error('MD.pixelSize_ is required.'); end
if isempty(MDtimeInterval_); error('MD.timeInterval_ is required.'); end
                                                    

disp(['== MDpixelSize_: ', num2str(MDpixelSize_), ' =='])
disp(['== MDtimeInterval_: ', num2str(MDtimeInterval_), ' =='])

%%  .txt (export comma delimited files)

fname0 = ['Chan', num2str(iChan)];

dlmwrite(fullfile(figuresDir, [fname0, '_velmap_outl.txt']), velmap_outl, 'precision', 8)
dlmwrite(fullfile(figuresDir, [fname0, 'imvelocitymap.txt']), imvelocitymap, 'precision', 8)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plots 

%%  Raw non-smoothActivityMap prot/act maps
%gaussianFilter = fspecial('gaussian', 13, 3);

%smParam = 1

inputmap = velmap;

%filteredmap = smoothActivityMap(velmap, 'SmoothParam', smParam, 'UpSample', 1);
fvelraw = figure('Visible', figFlag);  
figtmp = imagesc(inputmap);
title(chan0Title)
colorbar;colormap(jet) 

figtmp.AlphaData = 1-isnan(inputmap);
axis xy;xlabel('Time (s)');ylabel('Window')
ax = gca;
curTick = ax.XTick;
ax.XTickMode = 'manual';
ax.XTick = curTick+1;
ax.XTickLabel = (curTick)*MDtimeInterval_;

%%
saveas(fvelraw, fullfile(figuresDir, ['/raw', fname0, 'Map.png']), 'png')


%%  outl non-smoothActivityMap prot/act maps
%gaussianFilter = fspecial('gaussian', 13, 3);

%smParam = 1

inputmap = velmap_outl;

%filteredmap = smoothActivityMap(velmap, 'SmoothParam', smParam, 'UpSample', 1);
fvelraw = figure('Visible', figFlag);  
figtmp = imagesc(inputmap);
title(chan0Title)
colorbar;colormap(jet) 

figtmp.AlphaData = 1-isnan(inputmap);
axis xy;xlabel('Time (s)');ylabel('Window')
ax = gca;
curTick = ax.XTick;
ax.XTickMode = 'manual';
ax.XTick = curTick+1;
ax.XTickLabel = (curTick)*MDtimeInterval_;

%%
saveas(fvelraw, fullfile(figuresDir, ['outl_', fname0, 'Map.png']), 'png')



%%  Run:
%%  Boxplots after truncation and outlier detection

%
velBoxTime = figure('Visible', figFlag);
boxplot(velmap_outl, 'whisker', Inf);
title(chan0Title)
xlabel('time frame')
set(gca, 'XTick', 10:10:tmax)
set(gca, 'XTickLabel', 10:10:tmax)
 

velBoxWin = figure('Visible', figFlag);  
boxplot(velmap_outl', 'whisker', Inf);
title(chan0Title)
xlabel('Window')
set(gca, 'XTick', 10:10:wmax)
set(gca, 'XTickLabel', 10:10:wmax)

%%
saveas(velBoxTime, fullfile(figuresDir, [fname0, 'BoxTime.png']), 'png')
saveas(velBoxWin, fullfile(figuresDir, [fname0, 'BoxWin.png']), 'png')
 

%%  Histogram of velocity. see CV= sm.std/sm.mean*100

sm = summary(velmap_outl(:));
title2 = ['m=', sprintf('%0.2f', sm.mean), ' std=', sprintf('%0.2f', sm.std)];
velHist = figure('Visible', figFlag);
histogram(velmap_outl(:));
title({chan0Name, title2});

tmp = velmap_outl(:);

vpos = tmp(tmp >=0);
vneg = tmp(tmp < 0);

sm = summary(vpos);
title2 = ['m=', sprintf('%0.2f', sm.mean), ' std=', sprintf('%0.2f', sm.std)];
vposHist = figure('Visible', figFlag);  histogram(vpos);
title({['Positive ', chan0Name], title2})

sm = summary(vneg);
title2 = ['m=', sprintf('%0.2f', sm.mean), ' std=', sprintf('%0.2f', sm.std)];
vnegHist = figure('Visible', figFlag);  histogram(vneg);
title({['Negative ', chan0Name], title2})

%%
saveas(velHist, fullfile(figuresDir, [fname0, 'Hist.png']), 'png')
saveas(vposHist, fullfile(figuresDir, [fname0, 'PosHist.png']), 'png')
saveas(vnegHist, fullfile(figuresDir, [fname0, 'NegHist.png']), 'png')


%% topomap topographMD
if strcmp(p.topograph, 'on')

iWinProc = MD.getProcessIndex('WindowingProcess',1,0);
nBandMax_ = MD.processes_{iWinProc}.nBandMax_;
%
topovec = mean(velmap_outl, 2, 'omitnan');
topoMap = nan(wmax, nBandMax_);
topoMap(:, 1) = topovec;

title0 = chan0Title;

topomapFig = topographMD(MD, tmax, 1, topoMap, title0, figFlag);

%%
saveas(topomapFig, fullfile(figuresDir, [fname0, 'topomapFig.png']), 'png')

end

%%  smoothActivityMap prot/act maps

smParam = 0.8;

inputmap = velmap_outl;

filteredmap = smoothActivityMap(inputmap, 'SmoothParam', smParam, 'UpSample', 1);
fvel = figure('Visible', figFlag);
figtmp = imagesc(filteredmap);
title(chan0Title)
colorbar;colormap(jet)

figtmp.AlphaData = 1-isnan(inputmap);
axis xy;xlabel('Time (s)');ylabel('Window')
ax = gca;
curTick = ax.XTick;
ax.XTickMode = 'manual';
ax.XTick = curTick+1;
ax.XTickLabel = (curTick)*MDtimeInterval_;


%% 
saveas(fvel, fullfile(figuresDir, [fname0, 'Map.png']), 'png')

 
%%  Means plot

inputmap = velmap_outl;

% along time frame
cmeans = mean(inputmap, 1, 'omitnan');
cmeansZ = nanZscore(cmeans) ;

timeFr = 1:tmax;
timeAxis = (timeFr-1)*MDtimeInterval_;

meansTime = figure('Visible', figFlag);
plot(timeAxis, cmeansZ);
xlabel('Time (s)'); ylabel('Z-score')
title('Means plot')
legend(chan0Name, 'Location','northoutside','Orientation','horizontal')

% along windows
cmeans = mean(inputmap, 2, 'omitnan');
cmeansZ = nanZscore(cmeans) ;

winIndex = 1:wmax;

meansWin = figure('Visible', figFlag);
plot(winIndex, cmeansZ);
title('Means plot')
xlabel('Window'); ylabel('Z-score')
legend(chan0Name, 'Location','northoutside','Orientation','horizontal')


%% 
saveas(meansTime, fullfile(figuresDir, ['/meansTime', fname0, '.png']), 'png')
saveas(meansWin, fullfile(figuresDir, ['/meansWin', fname0, '.png']), 'png')


%%  TS plots for sampled 6 windows

inputmap = velmap_outl;
indNotAllNaN = find(~all(isnan(inputmap), 2));
ind0 = round(linspace(1, numel(indNotAllNaN), 6));
winInd = indNotAllNaN(ind0);

legend1 = {['win', num2str(winInd(1))], ['win', num2str(winInd(2))], ['win', num2str(winInd(3))]};

tsplots1 = figure('Visible', figFlag);
plot(timeAxis, inputmap(winInd(1:3), :))
xlabel('Time (s)'); ylabel(chan0Name)
title([chan0Title, ' example1'])
legend(legend1, 'Location','northoutside','Orientation','horizontal')

legend2 = {['win', num2str(winInd(4))], ['win', num2str(winInd(5))], ['win', num2str(winInd(6))]};

tsplots2 = figure('Visible', figFlag);
plot(timeAxis, inputmap(winInd(4:6), :))
xlabel('Time (s)'); ylabel(chan0Name)
title([chan0Title, ' example2'])
legend(legend2, 'Location','northoutside','Orientation','horizontal')

%% 
saveas(tsplots1, fullfile(figuresDir, ['tsplots1_', fname0, '.png']), 'png')
saveas(tsplots2, fullfile(figuresDir, ['tsplots2_', fname0, '.png']), 'png')


%%  spatial/temporal AutoCorr 1
% Select a map & change the figure name below
%map = imchan1map(40:126, 30:200);
% imvelocitymap

[acmap, corrMat, meanAutocorr] = ...
    TimeSpaceAutoCorPlot(imvelocitymap, chan0Name, MDtimeInterval_);

%%
saveas(acmap, fullfile(figuresDir, ['acmap', fname0, '.png']), 'png')
saveas(corrMat, fullfile(figuresDir, ['corrMat', fname0, '.png']), 'png')
saveas(meanAutocorr, fullfile(figuresDir, ['meanAutocorr', fname0, '.png']), 'png')


%% draw autoCorr curves (means) acCurve = autoCorrCurvePermTest

%numPerm = 1000;
%parpoolNum = 6;
%rseed = 'shuffle';

[acCurve, Avg_autocor] = autoCorrCurvePermTest(imvelocitymap, chan0Name, MDtimeInterval_, ...
                p.numPerm, p.parpoolNum, p.rseed);

%%
saveas(acCurve, fullfile(figuresDir, ['acCurve_', fname0, '.png']), 'png')
%%  03/23/2017
save(fullfile(figuresDir, [fname0, '_Avg_autocor_Vel.mat']), 'Avg_autocor')

%%  adftest map

if p.adf == 1
    mapForAdf = imvelocitymap(:, ~all(isnan(imvelocitymap)));
    
    [~, adfMap] = nanAdfTestMap(mapForAdf, chan0Name, 0.8);

%%
    saveas(adfMap, fullfile(figuresDir, ['adfMap', fname0, '.png']), 'png')  

end


%%  Coefficient of variation (sd/|mean|)

map = imvelocitymap;

% temporal
stdTemp = std(map, [], 2, 'omitnan');
%meanVec = mean(map, 2, 'omitnan');
%temporalCV = stdVec./abs(meanVec)

% spatial
stdSpat = std(map, [], 1, 'omitnan');

stdFull = [stdTemp ; stdSpat'];

tmax_ = size(map, 2); 
wmax_ = size(map, 1);

grChar = cell(wmax_+tmax_, 1);
for k = 1:wmax_ 
    grChar{k} = 'Temporal';
end
for k = 1:tmax_
    grChar{wmax_+k} = 'Spatial';
end

fvariation = figure('Visible', figFlag);
boxplot(stdFull, grChar)
title('Standard deviation (nm/sec)')

%%
saveas(fvariation, fullfile(figuresDir, ['CV', fname0, '.png']), 'png')


%%  checkWindowJump

winTrack = checkWindowJump(MD, 'figFlag', p.figFlag);

%
saveas(winTrack, fullfile(figuresDir, 'window1Trajectory.png'), 'png');


%%
disp('====End of mapDescriptives_Vel====')
    
end


