function MLsummary_XcorrCurvesVelAcf(ML, iChan1, iChan2, chan1Name, chan2Name, ...
    maxLayer, analNameAcf, analNameXcf, varargin)
% MLsummary_XcorrCurvesVelAcf Collect/Summarize cross correlation curves
% for movies computed and saved by mapXcorrCurvePermutation.m. 
% It computed 95% confidence bands of the Xcorrelations at each lag using 2*SEM 
% under the assumption that movies are independent replicates.
%
% Usage:
%       MLsummary_XcorrCurvesVelAcf(ML, 1, 0, 'Actin', 'Vel', 3, ...
%           'mapDescriptives', 'mapCrossCorr')
%
% Input:
%       ML          - a movieList object
%       iChan1      - the 1st channel index
%       chan1Name   - a short name for channel1.
%       iChan2      - the 2nd channel index
%       chan2Name   - a short name for channel2.
%       maxLayer    - maximum layer to be analyzed 
%       analNameAcf - the folder name for output from mapDescriptives_Vel.m
%                     to collect acf curves
%       analNameXcf - the folder name for output from
%                     mapXcorrCurvePermutation.m to collect xcf curves
%
% Output: .png/.fig/.mat files are saved in the ML.outputDirectory_/acfXcf_Ch1Ch0 by default. 
%       
% Option:
%       outDirName  - Specify a name of output directory.
%
% Jungsik Noh, 2017/05/17


% load ML, example MD
ML.getMovies()
md1 = ML.getMovie(1);
[~, cellLab0, ~] = fileparts(md1.outputDirectory_);

disp(['The label of the 1st movie will be the folder name for movieData.mat: ', cellLab0])

%%
MDs = ML.getMovies();
num = numel(MDs);

ch1ActmapName = chan1Name;
ch2ActmapName = chan2Name;

ch1fname = ['Ch', num2str(iChan1)];
ch2fname = ['Ch', num2str(iChan2)];
ch12fname = [ch1fname, ch2fname];

ip = inputParser; 
ip.addParameter('outDirName', ['acfXcf_', ch1fname, ch2fname]);
ip.addParameter('timeInterval', md1.timeInterval_);

parse(ip, varargin{:})
p = ip.Results;

fname_xcorrMat = ['xcorrMat_', ch1fname, ch2fname, '.mat'];
MDtimeIntvl = p.timeInterval;

%analNameAcf = 'mapDescriptives_0516_Imp0L4';  
%analNameXcf = 'mapCrossCorr_0510_Imp0L3';  
%% setting up parameters

%outDir = fullfile(ML.outputDirectory_, 'acfXcf_0515_Ch2Ch0');   % analName
outDir = fullfile(ML.outputDirectory_, p.outDirName);
if ~isdir(outDir); mkdir(outDir); end

%
acfvecsize = nan(num, 1);
for i = 1:num
    md = MDs{i};
    mdDir = md.outputDirectory_;
    %%%% input
    load(fullfile(mdDir, analNameAcf, 'Avg_autocor_Vel.mat')); % Avg_autocor
    acfvecsize(i) = size(Avg_autocor, 2);
end

lagSizeAcf = min(acfvecsize);
% or manually assign a value like: lagSizeAcf = 81;
disp(['Lag size (frames) of ACF: ', num2str(lagSizeAcf-1), ' frames'])

%
xcfmatsize = nan(num, 1);
for i = 1:num
    md = MDs{i};
    mdDir = md.outputDirectory_;
    %%%% input
    load(fullfile(mdDir, analNameXcf, fname_xcorrMat));  % xcorrMat
    xcfmatsize(i) = size(xcorrMat{1}, 2);
end
lagSizeXcf = min(xcfmatsize);
% or manually assign a value like: lagSizeXcf = 81;
disp(['Lag size (frames) of XCorrCurves: ', num2str(lagSizeXcf-1), ' frames'])



%%  Acf_vel

MDs = ML.getMovies();
num = numel(MDs);
%num = num    %%%  input
cellLabels = cell(num, 1);

acfArr = nan(num, lagSizeAcf);

%
for i = 1:num
    md = MDs{i};
    mdDir = md.outputDirectory_;
    [folderName, cellLab0, ~] = fileparts(mdDir);
    %cellName = [folderName((end-4):end), '_', cellLab0];   %% Adjust name extraction
    %cellName = [cellLab0(1:14)];   %% Adjust name extraction
    cellName = [cellLab0(1:end)];    % Generic cellName
    %[~, cellName] = fileparts(folderName);

    %%%% input
    load(fullfile(mdDir, analNameAcf, 'Avg_autocor_Vel.mat'));   

    %%%%
    
    cellLabels{i} = cellName;
    copyfile(fullfile(mdDir, analNameAcf, 'acCurveChan0.png'), ...
        fullfile(outDir, [cellLabels{i}, '_acCurveChan0.png']) )    
    
    %xcorrMat = Avg_autocor;
        %xcmean = mean(xcorrMat_tmp{indL}, 1, 'omitnan');
        xcmean = Avg_autocor;
        %tmplen = numel(xcmean);
        if (numel(xcmean) < lagSizeAcf)
            xcmeanMiddle = [xcmean, nan(1, lagSizeAcf - numel(xcmean))];
        else
            xcmeanMiddle = xcmean(1:lagSizeAcf);
        end
            
        acfArr(i, :) = reshape(xcmeanMiddle, 1, []);
    
end

totavgXcorrCurve = squeeze(mean(acfArr, 1, 'omitnan'));
 
 

%%
tLag = [0:lagSizeAcf-1].* MDtimeIntvl;

    fAcf = figure();    %'Visible', p.figFlag);   %%  name
    
    acmeanPerCellindL = squeeze(acfArr(:, :))';
    
    plot(tLag, acmeanPerCellindL, 'Color', [.5 .5 .5])
    h1 = refline([0,0]); h1.Color = [.5 .5 .5];
    %h = line([lagMax+1, lagMax+1], ylim); h.Color = [.5 .5 .5];
    hold on
    p1 = plot(tLag, totavgXcorrCurve, 'r');
    p1.LineWidth = 2;
    
    title('Vel')
    xlabel('Time lag (s)');ylabel('Correlation')
    set(gca, 'XGrid', 'on')

%%%%  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
ylim([-0.3, 0.3])


%% saveas

saveas(fAcf, fullfile(outDir, ['/acfVel', '.png']), 'png')
saveas(fAcf, fullfile(outDir, ['/acfVel', '.fig']), 'fig')





%%  Xcf

MDs = ML.getMovies();
num = numel(MDs);

%
num = num   %%%  input

cellLabels = cell(num, 1);

xcorrArr = nan(num, lagSizeXcf, maxLayer);
 
 
for i = 1:num
    md = MDs{i};
    mdDir = md.outputDirectory_;
    [folderName, cellLab0, ~] = fileparts(mdDir);
    %cellName = [folderName((end-4):end), '_', cellLab0];   %% Adjust name extraction
    cellName = [cellLab0(1:end)];    % Generic cellName
    %[~, cellName] = fileparts(folderName);

    %%%% input
    load(fullfile(mdDir, analNameXcf, fname_xcorrMat));   % xcorrMat_tmp
    %xcorrMat = xcorrMat_tmp;
    %%%%
    %pooledPermXcorrMean{i} = nan([size(permXcorrMean{1}), maxLayer]); 
    
    cellLabels{i} = cellName;
    copyfile(fullfile(mdDir, analNameXcf, ['xcCurve_', ch12fname, '_permT.png']), ...
        fullfile(outDir, [cellLabels{i}, '_xcorrCurve_', ch12fname, '_permT.png']) )    
    
    for indL = 1:maxLayer;
        %xcmean = mean(xcorrMat{indL}, 1, 'omitnan');
        xcmean = smoothingSplineCorMap(xcorrMat{indL});
        tmplen = numel(xcmean);
        xcmeanMiddle = xcmean(1+(tmplen - lagSizeXcf)/2:lagSizeXcf+(tmplen - lagSizeXcf)/2);
        xcorrArr(i, :, indL) = reshape(xcmeanMiddle, [], 1);
    end

end


totavgXcorrCurve = squeeze(mean(xcorrArr, 1, 'omitnan'));


%
save(fullfile(outDir, 'aggregatedXcorr.mat'), ...
      'cellLabels', 'xcorrArr', 'ch1ActmapName', ...
      'ch2ActmapName', 'maxLayer', 'analNameXcf', 'lagSizeXcf')



%%

lagMax = (lagSizeXcf - 1)/2;
lagGrid = floor(lagMax/2);
xcmapXtick = [1, 1+lagMax-lagGrid, 1+lagMax, 1+lagMax+lagGrid, 1+2*lagMax];
xcmapXticklabel = round( (xcmapXtick-1-lagMax)*MDtimeIntvl, 2);

%colOrd = get(gca,'ColorOrder');
%colOrd1 = colOrd(1:5, :);                   % Revise later together with legend.

f1 = cell(maxLayer,1);

for indL = 1:maxLayer

    f1{indL} = figure();    %'Visible', p.figFlag);   %%  name
    
    xcmeanPerCellindL = squeeze(xcorrArr(:, :, indL))';
    
    plot(1:lagSizeXcf, xcmeanPerCellindL, 'Color', [.5 .5 .5])
    h1 = refline([0,0]); h1.Color = [.5 .5 .5];
    %h = line([lagMax+1, lagMax+1], ylim); h.Color = [.5 .5 .5];
    hold on
    p1 = plot(1:lagSizeXcf, totavgXcorrCurve(:, indL), 'r');
    p1.LineWidth = 2;
    
    title(['xcorr(', ch1ActmapName, '_{t+Lag}, ', ch2ActmapName, '_t) -', num2str(indL), 'L'])
    xlabel('Time lag (s)');ylabel('Correlation')
%    legend('1L', '2L', '3L', '4L',  'Location','northoutside','Orientation','horizontal');
    set(gca, 'XTick', xcmapXtick)
    set(gca, 'XTickLabel', {xcmapXticklabel})
    %h = refline([0 0]); h.Color = [.5 .5 .5];
    %h = refline([Inf 0]); h.Color = [.5 .5 .5];
    set(gca, 'XGrid', 'on')
    
    hold on 

end


%% saveas
for indL=1:maxLayer
saveas(f1{indL}, fullfile(outDir, ['/xcorrSummary_', num2str(indL), 'L.png']), 'png')
saveas(f1{indL}, fullfile(outDir, ['/xcorrSummary_', num2str(indL), 'L.fig']), 'fig')
end


 
%%  integrate movies SEM

lagMax = (lagSizeXcf - 1)/2;
lagGrid = floor(lagMax/2);
xcmapXtick = [1, 1+lagMax-lagGrid, 1+lagMax, 1+lagMax+lagGrid, 1+2*lagMax];
xcmapXticklabel = round( (xcmapXtick-1-lagMax)*MDtimeIntvl, 2);

%colOrd = get(gca,'ColorOrder');
%colOrd1 = colOrd(1:5, :);                   % Revise later together with legend.

f2 = cell(maxLayer,1);

for indL = 1:maxLayer

    f2{indL} = figure();    %'Visible', p.figFlag);   %%  name
    
    xcMVecsIndL = squeeze(xcorrArr(:, :, indL));
    xcTotM = mean(xcMVecsIndL, 1, 'omitnan');
    sem = std(xcMVecsIndL, [], 1, 'omitnan')./sqrt(num);
    
    s1 = shadedErrorBar(1:lagSizeXcf, xcTotM, 2*sem, '-r', 1);
    s1.mainLine.LineWidth = 2;
    
    hold on
    h1 = refline([0,0]); h1.Color = [.5 .5 .5];
    
    title(['xcorr(', ch1ActmapName, '_{t+Lag}, ', ch2ActmapName, '_t) -', num2str(indL), 'L'])
    xlabel('Time lag (s)');ylabel('Correlation')
%    legend('1L', '2L', '3L', '4L',  'Location','northoutside','Orientation','horizontal');
    set(gca, 'XTick', xcmapXtick)
    set(gca, 'XTickLabel', {xcmapXticklabel})
    %h = refline([0 0]); h.Color = [.5 .5 .5];
    %h = refline([Inf 0]); h.Color = [.5 .5 .5];
    set(gca, 'XGrid', 'on')
    
    hold on 

end


%% saveas
for indL=1:maxLayer
saveas(f2{indL}, fullfile(outDir, ['/xcorrMean_', num2str(indL), 'L.png']), 'png')
saveas(f2{indL}, fullfile(outDir, ['/xcorrMean_', num2str(indL), 'L.fig']), 'fig')
end




%%  CMLags

MDs = ML.getMovies();
num = numel(MDs);


%%

cellLabels = cell(num, 1);
CorMaxLags = nan(num, 2*maxLayer);
%rownames = {'cm1L', 'lag1L', 'cm2L', 'lag2L', 'cm3L', 'lag3L'};
%rownames = {'CM1L', 'lag1L'};
rownames = cell(1, 2*maxLayer);
for indL = 1:maxLayer
    rownames{1 + 2*(indL-1)} = ['CM', num2str(indL), 'L'];
    rownames{2*indL} = ['lag', num2str(indL), 'L'];
end


for i = 1:num
    md = MDs{i};
    mdDir = md.outputDirectory_;
    [folderName, cellLab0, ~] = fileparts(mdDir);
    %cellName = [folderName((end-4):end), '_', cellLab0];   %% Adjust name extraction
    cellName = [cellLab0(1:end)];    % Generic cellName
    %[~, cellName] = fileparts(folderName);

    %%%% input
    %analName = 'mapCrossCorr_0207_Imp1L3';
    load(fullfile(mdDir, analNameXcf, fname_xcorrMat));
    %%%%
    
    cmCell = xcorrMatToCMLag(xcorrMat);
    cmMat = cell2mat(cmCell)';
    cmVec = reshape(cmMat, 1, []);
    
    for j = 1:numel(cmVec)
        CorMaxLags(i, j) = cmVec(j);
    end
    cellLabels{i} = cellName;
end

T0 = table(cellLabels);
T1 = array2table(CorMaxLags);
T1.Properties.VariableNames = rownames;
T = [T0, T1];

%%
save(fullfile(outDir, 'CorrMaximaLag_table.mat'), 'T', 'cellLabels', 'CorMaxLags')



%% cor maxima/minima plot

fcm = cell(maxLayer, 1);

ind0 = 1:maxLayer;
cors = CorMaxLags(:, 1+2*(ind0-1));
lags = CorMaxLags(:, 2*ind0);

corsMin = min(min(cors(:)), -0.01) * 1.2;
corsMax = max(max(cors(:)), 0.01) * 1.2;
lagsMin = min(lags(:)) - 3;
lagsMax = max(lags(:)) + 3;


for indL=1:maxLayer

    r1 = CorMaxLags(:, 1+2*(indL-1));

    l1 = CorMaxLags(:, 2*indL) .*MDtimeIntvl;

    % xlim0 = [-15, 50]
    % ylim([-0.1 0.15]) 

    fcm{indL} = figure;
    scatter(l1, r1, [], 'r', 'filled')
    text(l1+0.5, r1+rand(num,1)*0.02-0.01, cellLabels)
    ylabel('Corr Maximum/Minmum')
    xlabel('Lag (relative to vel) (sec)')
    
    rMean = mean(r1);
    lMean = mean(l1);
    lMedian = median(l1);
    rMedian = median(r1);
    title1 = ['Layer ', num2str(indL)];
    title2 = ['avgLag: ', num2str(lMean), '(', num2str(lMedian), ')', ' avgCor: ', num2str(rMean), '(', num2str(rMedian), ')'];
    title({title1, title2})
    hold on
    scatter(lMean, rMean, [], '+')

    ylim([corsMin, corsMax])
    xlim0 = [lagsMin lagsMax].* MDtimeIntvl;
    xlim(xlim0)

    h=refline([0,0]); %h.Color = [.5 .5 .5];
    h = line([0, 0], ylim); %h.Color = [.5 .5 .5];

end



%% saveas
for indL=1:maxLayer
saveas(fcm{indL}, fullfile(outDir, ['/CMLag', num2str(indL), 'L.png']), 'png')
saveas(fcm{indL}, fullfile(outDir, ['/CMLag', num2str(indL), 'L.fig']), 'fig')
end



%%
disp('==== MLsummary_XcorrCurvesVelAcf is finished!! ====')
disp('==:) for i=1:40; close(figure(i)); end')

end
