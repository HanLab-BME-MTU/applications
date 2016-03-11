function analyzeWindowSpatialVariation(ML)
%UNDER CONSTRUCTION, BETA, BLA BLA BLA
%RIGHT NOW THIS IS ONLY WORKING ON SINGLE_FRAME_MOVIES!!!!!!!! STAY AWAY!!!

%% ---------------- Parameters ----------------- %%

makePlots = true;
movOutName = 'windowing spatial variation analysis';
pOpt = {'-r300',...% dpi = 300
        '-depsc2'};% use eps format
pOptTIFF = {'-r100','-dtiff'};%150 dpi for TIF format since this is usally just used for emailing to Bob!

cMapSat = 3;%Percentage saturation on colormaps

%% -------------- Pre-Processing --------------- %%

%outDir = [ML.outputDirectory_ filesep ML.movieListFileName_(1:end-4)];
outDir = ML.outputDirectory_;


if ~exist(outDir,'dir')
    mkdir(outDir);
end

nMov = numel(ML.movies_);
sampAllMov = cell(nMov,1);
nChan = zeros(nMov,1);
sampSize = nan(nMov,2);
hasFeat = false(nMov,1);
allWin = cell(nMov,1);
iWinSampProc = zeros(nMov,1);
iWinProc = zeros(nMov,1);
movOutDirs = cell(nMov,1);

firstFound = true;


for j = 1:nMov
    
    
    
    %Set up the output for this movie
    movOutDirs{j} = [ML.movies_{j}.outputDirectory_ filesep movOutName];
    mkClrDir(movOutDirs{j})   

    %TEMP - right now we know they're all 4 chans!
    nChan(j) = numel(ML.movies_{j}.channels_);
    sampAllMov{j} = cell(nChan(j),1);
    
    %Load the samples for this movie
    iWinSampProc(j) = ML.movies_{j}.getProcessIndex('WindowSamplingProcess',1,0);
    for k = 1:nChan(j)
        sampAllMov{j}{k} = ML.movies_{j}.processes_{iWinSampProc(j)}.loadChannelOutput(k);
    end           
            
    %TEMP - get from process, initialize, don't use cell etc etc. 
    sampSize(j,:) = size(sampAllMov{j}{k}.avg);   
    
    %Load the feature identification if available
    featFile = [ML.movies_{j}.outputDirectory_ filesep  'manual morphological feature identification.mat'];
    if exist(featFile,'file')
        hasFeat(j) = true;
        tmp = load(featFile,'featIdent');
        featIdent(j) = tmp.featIdent;
        
        if firstFound
            %Get info about morpho features
            featNames = fieldnames(featIdent);
            featNames = featNames(~strcmp('Time',featNames));
            nFeat = numel(featNames);
            iBackFeat = find(strcmp('Back',featNames));%TEMP??        
            featCloseWin = zeros(nMov,nFeat);
            firstFound = false;
        end       
    end    
    
    %And the windows - we need these for figures, feature matching etc
    iWinProc(j) = ML.movies_{j}.getProcessIndex('WindowingProcess',1,0);
    allWin{j} = ML.movies_{j}.processes_{iWinProc(j)}.loadChannelOutput(1);
    
    featCloseWin(j,iBackFeat) = findClosestWindow(allWin{j},featIdent(j).(featNames{iBackFeat}));
    
    %Permute the matrices for this channel so that the origin is at the
    %back for sample alignment
    iPerm{j} = [featCloseWin(j,iBackFeat):sampSize(j,1) 1:(featCloseWin(j,iBackFeat)-1)];%Index for permutation
    allWin{j} = allWin{j}(iPerm{j});%Permute the windows to match

%     iP = 1:sampSize(j,1);%Get indices for permuting the other indices. Is there an easier way to do this?
%     iPP = iP(iPerm{j}(end:-1:1));
%     iLam{j} = iPP(iLam{j});
%     iOther{j} = iPP(iOther{j});
    
    %Get the closest window for each morphological feature
    for k = 1:nFeat
        featCloseWin(j,k) = findClosestWindow(allWin{j},featIdent(j).(featNames{k}));
    end
    
    %Get indices for lamellar and "other" (~back) regions TEMP - THIS IS NOT GENERIC!!WHO
    %GIVES A SHIT THOUGH???
    if featCloseWin(j,2) < featCloseWin(j,3)%If the lamellar region crosses the origin
        iLam{j} = [featCloseWin(j,3):sampSize(j,1) 1:featCloseWin(j,2)];
        iOther{j} = featCloseWin(j,2)+1:featCloseWin(j,3)-1;                
    else%If the "other" region crosses the origin
        iLam{j} = featCloseWin(j,3):featCloseWin(j,2);
        iOther{j} = [featCloseWin(j,2)+1:sampSize(j,1) 1:featCloseWin(j,3)-1];
    end
    lamSize(j) = numel(iLam{j});
    othSize(j) = numel(iOther{j});
    
    
    
    
    if makePlots
        %Make Region Labelling / Permutation Plot
        %currFig = fsFigure(.33);
        currFig = imageViewer(ML.movies_{j},'ChannelIndex',3,'Saturate',.03);
        hold on
        plotWindows(allWin{j},{'k','FaceAlpha',0},5)        
        plotWindows(allWin{j}(iLam{j}),{'b','FaceAlpha',.2},0)
        plotWindows(allWin{j}(iOther{j}),{'r','FaceAlpha',.2},0)    
        text(featIdent(j).Back(1),featIdent(j).Back(2),'Back','Color','y','FontSize',16)
        %plot(featIdent(j).Back(1),featIdent(j).Back(2),'rx')
        text(featIdent(j).Lam1(1),featIdent(j).Lam1(2),'Lamella 1','Color','y','FontSize',16)
        text(featIdent(j).Lam2(1),featIdent(j).Lam2(2),'Lamella 2','Color','y','FontSize',16)
        hgsave(currFig,[movOutDirs{j} filesep 'region labelling figure.fig'])
        print(currFig,[movOutDirs{j} filesep 'region labelling figure.tif'],pOptTIFF{:})
    end    
    
    
end

minSize = min(sampSize,[],1);
minLamSize = min(lamSize);
minOthSize = min(othSize);


sampAllIntBothAvg = nan([minSize nMov max(nChan)]);
sampAllIntParaAvg = nan([minSize nMov max(nChan)]);
% sampAllIntPerpAvg = nan([maxSize nMov nChan]);
% sampAllIntParaAvg = nan([maxSize nMov nChan]);

sampLamIntBothAvg = nan([minLamSize minSize(2) nMov max(nChan)]);
sampLamIntParaAvg = nan([minLamSize minSize(2) nMov max(nChan)]);
sampOthIntBothAvg = nan([minOthSize minSize(2) nMov max(nChan)]);
sampOthIntParaAvg = nan([minOthSize minSize(2) nMov max(nChan)]);

%Go through each, align samples and interpolate to match the minimum size

for j = 1:nMov
               
    for k = 1:nChan(j)        
                       
        %Permute the samples for each channel also
        sampAllMov{j}{k}.avg = sampAllMov{j}{k}.avg(iPerm{j},:);
                            
        [Xi,Yi] = meshgrid(linspace(1,sampSize(j,2),minSize(2)),linspace(1,sampSize(j,1),minSize(1)));                       
        sampAllIntBothAvg(:,:,j,k) = interp2(sampAllMov{j}{k}.avg,Xi,Yi);                        
        [Xi,Yi] = meshgrid(1:sampSize(j,2),linspace(1,sampSize(j,1),minSize(1)));                       
        sampAllIntParaAvg(:,1:sampSize(j,2),j,k) = interp2(sampAllMov{j}{k}.avg,Xi,Yi);       
        
        
        [Xi,Yi] = meshgrid(linspace(1,sampSize(j,2),minSize(2)),linspace(1,lamSize(j),minLamSize));                        
        sampLamIntBothAvg(:,:,j,k) = interp2(sampAllMov{j}{k}.avg(iLam{j},:),Xi,Yi);                        
        [Xi,Yi] = meshgrid(1:sampSize(j,2),linspace(1,lamSize(j),minLamSize));                        
        sampLamIntParaAvg(:,1:sampSize(j,2),j,k) = interp2(sampAllMov{j}{k}.avg(iLam{j},:),Xi,Yi);       
        
        [Xi,Yi] = meshgrid(linspace(1,sampSize(j,2),minSize(2)),linspace(1,othSize(j),minOthSize));                        
        sampOthIntBothAvg(:,:,j,k) = interp2(sampAllMov{j}{k}.avg(iOther{j},:),Xi,Yi);
        [Xi,Yi] = meshgrid(1:sampSize(j,2),linspace(1,othSize(j),minOthSize));                        
        sampOthIntParaAvg(:,1:sampSize(j,2),j,k) = interp2(sampAllMov{j}{k}.avg,Xi,Yi);       
                        
                    
    end               
        
end





%% ------------- All-Cell Analysis ---------------- %%

scanFigPara = figure;
hold on
scanFigPerp = figure;
hold on
scanFigPerpParaInterp = figure;
hold on

scanFigParaLam = figure;
hold on
scanFigPerpLam = figure;
hold on
scanFigPerpLamParaInterp = figure;
hold on

scanFigParaOth = figure;
hold on
scanFigPerpOth = figure;
hold on
scanFigPerpOthParaInterp = figure;
hold on

chanCols = jet(max(nChan));
chanCols = chanCols(:,end:-1:1);

%TEMP!
chanNames = {'p-Erk','Wave2','Phalloidin','DAPI'};

for k = 1:max(nChan)
    
    % --- Whole-Cell ---- %
    
    wholeCellMeanMap = squeeze(nanmean(sampAllIntBothAvg(:,:,:,k),3));
    
    figure    
    imagesc(squeeze(nanmean(sampAllIntBothAvg(:,:,:,k),3))),axis image
    saturateImageColormap(gcf,cMapSat);
    ca = caxis;
    title([chanNames{k} '  Whole-Cell Interp Both'])
    hgsave([outDir filesep 'full spatial profile all cell average interp both directions channel ' num2str(k)])
    
    figure    
    imagesc(squeeze(nanmean(sampAllIntParaAvg(:,:,:,k),3))),axis image
    saturateImageColormap(gcf,cMapSat);
    ca = caxis;
    title([chanNames{k} '  Whole-Cell Interp Para'])
    hgsave([outDir filesep 'full spatial profile all cell average interp para only channel ' num2str(k)])
    
    
%     figure
%     polarplot3d(wholeCellMeanMap(:,round(end/2):-1:1)');view(2)    
%     caxis(ca);
%     title([chanNames{k} '  Whole-Cell'])
%     hgsave([outDir filesep 'full spatial profile all cell average interp both directions channel ' num2str(k)])
%     
    
    
    figure(scanFigPara)
    plot(squeeze(nanmean(nanmean(sampAllIntBothAvg(:,:,:,k),2),3)),'color',chanCols(k,:))
        
    figure(scanFigPerp)
    plot(squeeze(nanmean(nanmean(sampAllIntBothAvg(:,:,:,k),1),3)),'color',chanCols(k,:))
    figure(scanFigPerpParaInterp)
    plot(squeeze(nanmean(nanmean(sampAllIntParaAvg(:,:,:,k),1),3)),'color',chanCols(k,:))
    
    % --- Lamellar --- %
    
    figure    
    imagesc(squeeze(nanmean(sampLamIntBothAvg(:,:,:,k),3))),axis image
    saturateImageColormap(gcf,cMapSat);
    title([chanNames{k} '  lamellar areas Interp Both'])
    hgsave([outDir filesep 'lamellar spatial profile all cell average interp both directions channel ' num2str(k)])
    
    figure    
    imagesc(squeeze(nanmean(sampLamIntParaAvg(:,:,:,k),3))),axis image
    saturateImageColormap(gcf,cMapSat);
    title([chanNames{k} '  lamellar areas, Interp Para'])
    hgsave([outDir filesep 'lamellar spatial profile all cell average interp para only channel ' num2str(k)])
    
    
    figure(scanFigParaLam)
    plot(squeeze(nanmean(nanmean(sampLamIntBothAvg(:,:,:,k),2),3)),'color',chanCols(k,:))
        
    figure(scanFigPerpLam)
    plot(squeeze(nanmean(nanmean(sampLamIntBothAvg(:,:,:,k),1),3)),'color',chanCols(k,:))
    figure(scanFigPerpLamParaInterp)
    plot(squeeze(nanmean(nanmean(sampLamIntParaAvg(:,:,:,k),1),3)),'color',chanCols(k,:))
    
    % ---- Other ---- %
    
    figure    
    imagesc(squeeze(nanmean(sampOthIntBothAvg(:,:,:,k),3))),axis image
    saturateImageColormap(gcf,cMapSat);
    title([chanNames{k} '  non-lamellar areas, interp both'])
    hgsave([outDir filesep 'non lamellar spatial profile all cell average interp both directions channel ' num2str(k)])        
    
    figure    
    imagesc(squeeze(nanmean(sampOthIntParaAvg(:,:,:,k),3))),axis image
    saturateImageColormap(gcf,cMapSat);
    title([chanNames{k} '  non-lamellar areas, interp para'])
    hgsave([outDir filesep 'non lamellar spatial profile all cell average interp para only channel ' num2str(k)])        
    
    
    figure(scanFigParaOth)
    plot(squeeze(nanmean(nanmean(sampOthIntBothAvg(:,:,:,k),2),3)),'color',chanCols(k,:))
        
    figure(scanFigPerpOth)
    plot(squeeze(nanmean(nanmean(sampOthIntBothAvg(:,:,:,k),1),3)),'color',chanCols(k,:))        
    figure(scanFigPerpOthParaInterp)
    plot(squeeze(nanmean(nanmean(sampOthIntParaAvg(:,:,:,k),1),3)),'color',chanCols(k,:))        
    
end

figure(scanFigPara)
title('Parallel Scan, Whole-Cell')
legend(chanNames)
figure(scanFigPerp)
title('Perpindicular Scan, Whole-Cell')
legend(chanNames)
figure(scanFigPerpParaInterp)
title('Perpindicular Scan, Whole-Cell, Para Interp Only')
legend(chanNames)


figure(scanFigParaLam)
title('Parallel Scan, Lamellar Areas')
legend(chanNames)
figure(scanFigPerpLam)
title('Perpindicular Scan, Lamellar Areas')
legend(chanNames)
figure(scanFigPerpLamParaInterp)
title('Perpindicular Scan, Lamellar Areas, Para Interp Only')
legend(chanNames)

figure(scanFigParaOth)
title('Parallel Scan, Non-Lamellar Areas')
legend(chanNames)
figure(scanFigPerpOth)
title('Perpindicular Scan, Non-Lamellar Areas')
legend(chanNames)
figure(scanFigPerpOthParaInterp)
title('Perpindicular Scan, Non-Lamellar Areas, Para Interp Only')
legend(chanNames)


%Do the bands second so the legends work
% for k =  1:max(nChan)
%     
%     figure(scanFigPara)
%     plotTransparent(1:maxSize(1),squeeze(nanmean(nanmean(sampAllIntBothAvg(:,:,:,k),2),3)),...
%         squeeze(nanstd(nanmean(sampAllIntBothAvg(:,:,:,k),2),[],3)),chanCols(k,:),.2,0)
%     
%     figure(scanFigPerp)    
%     plotTransparent(1:maxSize(2),squeeze(nanmean(nanmean(sampAllIntBothAvg(:,:,:,k),1),3)),...
%         squeeze(nanstd(nanmean(sampAllIntBothAvg(:,:,:,k),1),[],3)),chanCols(k,:),.2,0)    
%     figure(scanFigPerpParaInterp)    
%     plotTransparent(1:maxSize(2),squeeze(nanmean(nanmean(sampAllIntParaAvg(:,:,:,k),1),3)),...
%         squeeze(nanstd(nanmean(sampAllIntParaAvg(:,:,:,k),1),[],3)),chanCols(k,:),.2,0)    
% 
%     figure(scanFigParaLam)
%     plotTransparent(1:maxLamSize(1),squeeze(nanmean(nanmean(sampLamIntBothAvg(:,:,:,k),2),3)),...
%         squeeze(nanstd(nanmean(sampLamIntBothAvg(:,:,:,k),2),[],3)),chanCols(k,:),.2,0)
%     
%     figure(scanFigPerpLam)    
%     plotTransparent(1:maxSize(2),squeeze(nanmean(nanmean(sampLamIntBothAvg(:,:,:,k),1),3)),...
%         squeeze(nanstd(nanmean(sampLamIntBothAvg(:,:,:,k),1),[],3)),chanCols(k,:),.2,0)    
%     figure(scanFigPerpLamParaInterp)    
%     plotTransparent(1:maxSize(2),squeeze(nanmean(nanmean(sampLamIntParaAvg(:,:,:,k),1),3)),...
%         squeeze(nanstd(nanmean(sampLamIntParaAvg(:,:,:,k),1),[],3)),chanCols(k,:),.2,0)    
%     
%     figure(scanFigParaOth)
%     plotTransparent(1:maxOthSize(1),squeeze(nanmean(nanmean(sampOthIntBothAvg(:,:,:,k),2),3)),...
%         squeeze(nanstd(nanmean(sampOthIntBothAvg(:,:,:,k),2),[],3)),chanCols(k,:),.2,0)
%     
%     figure(scanFigPerpOth)    
%     plotTransparent(1:maxSize(2),squeeze(nanmean(nanmean(sampOthIntBothAvg(:,:,:,k),1),3)),...
%         squeeze(nanstd(nanmean(sampOthIntBothAvg(:,:,:,k),1),[],3)),chanCols(k,:),.2,0)    
%     figure(scanFigPerpOthParaInterp)    
%     plotTransparent(1:maxSize(2),squeeze(nanmean(nanmean(sampOthIntParaAvg(:,:,:,k),1),3)),...
%         squeeze(nanstd(nanmean(sampOthIntParaAvg(:,:,:,k),1),[],3)),chanCols(k,:),.2,0)    
% 
% end

hgsave(scanFigPara,[outDir filesep 'along edge profile all cell average interp both directions all channels']);
hgsave(scanFigPerp,[outDir filesep 'away from edge profile all cell average interp both directions all channels']);
hgsave(scanFigPerpParaInterp,[outDir filesep 'away from edge profile all cell average interp para only all channels']);

hgsave(scanFigParaLam,[outDir filesep 'along edge profile lamellar average interp both directions all channels']);
hgsave(scanFigPerpLam,[outDir filesep 'away from edge profile lamellar average interp both directions all channels']);
hgsave(scanFigPerpLamParaInterp,[outDir filesep 'away from edge profile lamellar average interp para only all channels']);

hgsave(scanFigParaOth,[outDir filesep 'along edge profile oth average interp both directions all channels']);
hgsave(scanFigPerpOth,[outDir filesep 'away from edge profile oth average interp both directions all channels']);
hgsave(scanFigPerpOthParaInterp,[outDir filesep 'away from edge profile oth average interp para only all channels']);



figure
hold on
title('Coverage, Interp Both')
imagesc(squeeze(sum(~isnan(sampAllIntBothAvg(:,:,:,1)),3))),axis image,colorbar
hgsave([outDir filesep 'coverage all cells interp both directions'])

save([outDir filesep 'combined sample analysis.mat'], 'sampAllIntBothAvg','sampAllIntParaAvg','sampOthIntBothAvg','sampLamIntBothAvg',...
       'ML','allWin','chanNames','iLam','iOther','iPerm','sampSize','sampAllMov','wholeCellMeanMap');


