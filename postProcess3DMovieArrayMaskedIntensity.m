function postProcess3DMovieArrayMaskedIntensity(MA,p)

iProcChan = 1;

nMov = numel(MA);

pOpt = {'-r300',...% dpi = 300
        '-depsc2'};% use eps format
pOptTif = {'-r100',...% dpi = 100
        '-dtiff'};% use tiff format

outDir = p.OutputDirectory;    

useMicrons = true;

mkClrDir(p.OutputDirectory)

intAn = cell(nMov,1);
nFramesPerMov = nan(nMov,1);
nChanPer = nan(nMov,1);
pixelSizePer = nan(nMov,1);

if ~isfield(p,'ChannelIndex')
    p.ChannelIndex = [];
end
if ~isfield(p,'NormalizeIntensities')
    p.NormalizeIntensities = true;
end


%If we want to restrict analysis to first n frame to avoid differential
%bleaching effects and give fair comparision with fixed cell data
if ~isfield(p,'UseFrames')
    p.UseFrames = [];
end

curvConv = cell(nMov,1);
[intTypes intNames] = getIntTypeFields;
nIntTypes = numel(intTypes);

iIntUse = [1 4 5 8];%Don't show all the intensities to keep things simple

%Get limits and indices for curvature categories.
[curvCatKLims,curvCatHLims,curvCatNames,curvCatColors] = getCurvCategories;

nCurvCat = numel(curvCatNames);



%% ------------- Per-Movie Loading and Processing --------------- %%


for iMov = 1:nMov
    
    if isempty(MA(iMov).eventTimes_)
        nFramesPerMov(iMov) = MA(iMov).nFrames_;
    else
        nFramesPerMov(iMov) = MA(iMov).eventTimes_;
    end    
    
    if ~isempty(p.UseFrames)
        nFramesPerMov(iMov) = min(nFramesPerMov(iMov),p.UseFrames);
    end
    
    pixelSizePer(iMov) = MA(iMov).pixelSize_;
    
    if isempty(p.ChannelIndex)
        nChanPer(iMov) = numel(MA(iMov).channels_);
    else
        %If the user specifies we use only the first n channels
        nChanPer(iMov) = min(numel(MA(iMov).channels_),numel(p.ChannelIndex));
    end
    
    %So we can later check that the channels/fluorophores are in agreement
    chanPaths{iMov} = arrayfun(@(x)(x.channelPath_),MA(iMov).channels_,'Unif',0);
    
    iMAProc = MA(iMov).getProcessIndex('MaskedIntensity3DProcess',1 ,0);        
    assert(~isempty(iMAProc),['Movie ' num2str(iMov) ' has no intensity analysis!']);    
        
    intAn{iMov} = MA(iMov).processes_{iMAProc}.loadChannelOutput(iProcChan);
    
    %Remove any bad frames from the branch profiles
    intAn{iMov}.branchProfiles = intAn{iMov}.branchProfiles(1:nFramesPerMov(iMov));
    
    [curvTypes,curvNames,curvUnits,curvConv{iMov}] = getCurveTypeFields(pixelSizePer(iMov),useMicrons);
    nCurvTypes = numel(curvTypes);               
        
    if nFramesPerMov(iMov) > 1
    
        %TEEEEMPPP ADD MULTI-CHANNEL SUPPORT??? Probably wont see
        %multi-channel timelaps for years....
        ccPerMov(iMov,:,:) = intAn{iMov}.ccIntCurv(1,:,:);
        tIntPer(iMov) = MA(iMov).timeInterval_;                
    else
        tIntPer(iMov) = 0;
    end
        
    
    %Extract and combine per-frame data
                
    
    nSurfPtsPerFrame{iMov} = cellfun(@numel,{intAn{iMov}.branchProfiles(:).surfPixInd});
    nSurfPtsTot(iMov) = sum(nSurfPtsPerFrame{iMov});
    allCurvPerMov{iMov} = nan(nSurfPtsTot(iMov),nCurvTypes);    
    for j = 1:nCurvTypes
        allCurvPerMov{iMov}(:,j) = vertcat( intAn{iMov}.branchProfiles(:).(curvTypes{j}) ) .* curvConv{iMov}(j);            
    end
    
    allIntPerMov{iMov} = nan(nSurfPtsTot(iMov),nIntTypes,nChanPer(iMov));    
    for j = 1:nIntTypes    
        if ~isempty(p.ChannelIndex)
            tmp = vertcat(intAn{iMov}.branchProfiles(:).(intTypes{j}));     
            allIntPerMov{iMov}(:,j,:) = tmp(:,p.ChannelIndex);
        else
            allIntPerMov{iMov}(:,j,:) = vertcat(intAn{iMov}.branchProfiles(:).(intTypes{j}));     
        end
    end
    
    for iFrame = 1:nFramesPerMov(iMov)
        
        nBranchPerFrame{iMov}(iFrame) = numel(intAn{iMov}.branchProfiles(iFrame).branchTipPixelLenVsDepth);            
        hasLvD = cellfun(@(x)(~isempty(x) && ~any(isnan(x(:)))),intAn{iMov}.branchProfiles(iFrame).branchTipPixelLenVsDepth);                        
        
        avgBranchRadiusPerFrame{iMov,iFrame} = cellfun(@(x)(mean(x(:,2))),intAn{iMov}.branchProfiles(iFrame).branchTipPixelLenVsDepth(hasLvD)) .* MA(iMov).pixelSize_;                                
        for k = 1:nChanPer(iMov);
            
            avgBranchIntPerFrame{iMov,iFrame}(:,k) = cellfun(@(x)(mean(x(:,k))),intAn{iMov}.branchProfiles(iFrame).branchTipPixelInt(hasLvD));
            maxBranchIntPerFrame{iMov,iFrame}(:,k) = cellfun(@(x)(max(x(:,k))),intAn{iMov}.branchProfiles(iFrame).branchTipPixelInt(hasLvD));
            medBranchIntPerFrame{iMov,iFrame}(:,k) = cellfun(@(x)(median(x(:,k))),intAn{iMov}.branchProfiles(iFrame).branchTipPixelInt(hasLvD));

        end                
        branchLenPerFrame{iMov,iFrame} = cellfun(@(x)(max(x(:,1))),intAn{iMov}.branchProfiles(iFrame).branchTipPixelLenVsDepth(hasLvD)) .* MA(iMov).pixelSize_;
                        
        %Get sampling statistics for estimating sample overlap
        nFacesSampledPerFame{iMov}(iFrame) = sum(intAn{iMov}.branchProfiles(iFrame).nFacesPerCurvSamp);        
        currMG = MA(iMov).processes_{MA(iMov).getProcessIndex('MaskGeometry3DProcess',1,0)}.loadChannelOutput(1,iFrame);
        nFacesPerFrame{iMov}(iFrame) = size(currMG.SmoothedSurface.faces,1);
        %Now get the average number of times each face was sampled
        nResampFacesPerFrame{iMov}(iFrame) = nFacesSampledPerFame{iMov}(iFrame) / nFacesPerFrame{iMov}(iFrame);
        
    end
      
end

assert(numel(unique(nChanPer))==1,'All movies must have the same number of channels!!')
nChan = nChanPer(1);
if isempty(p.ChannelIndex)
    p.ChannelIndex = 1:nChan;
end

if numel(unique(tIntPer)) > 1 || numel(unique(nFramesPerMov)) > 1
    %Obviously this isn't strictly true, we're just using this constraint
    %for now.
    warning('All movies should have the same time-interval and number of frames to run combined cross-corr!!')    
end


%% ------------- Cross-Correlation Analysis  --------------- %%

tInt = tIntPer(1);    
nFrames = nFramesPerMov(1);
isTLapse = tInt ~= 0;

if isTLapse

    maxDelay = (size(ccPerMov,3)-1)/2;
    tLags = -tInt*maxDelay:tInt:tInt*maxDelay;

    for j = 1:nCurvTypes

        figure
        % waterfall(tLags,1:nMov,ccPerMov)
        % zlabel('Correlation')
        % ylabel('Time Delay, Seconds')
        % xlabel('Movie #')
        hold on
        for k = 1:nMov
            plot(tLags,squeeze(ccPerMov(k,j,:)),'color',rand(1,3))    
        end
        xlabel('Time Delay, Seconds')
        ylabel('Correlation')
        plot(xlim,[0 0],'--k')
        plot([ 0 0],ylim,'--k')
        title({['Cross-Correlation, Mean Cortical Intensity and ' curvTypes{j} ],...
                'Positive Delay means curvature follows intensity'})
        saveThatShit(['cross corr individual movies ' curvTypes{j}],p.OutputDirectory)



        [ccMean,ccCI] = correlationBootstrap(squeeze(ccPerMov(:,j,:))',ones(1,nMov)*1.96/sqrt(nFrames));


        figure
        plot(tLags,ccMean,'.-')
        hold on
        plot(tLags,ccCI(1,:),'--k')
        legend('Mean','95% C.I.')    
        plot(tLags,ccCI(2,:),'--k')
        plot(xlim,[0 0 ],'--k')
        plot([ 0  0],ylim,'--k')
        xlabel('Time Delay, Seconds')
        ylabel('Correlation')
        title({['Cross-Correlation, Mean Cortical Intensity and ' curvTypes{j} ],...
                'Positive Delay means curvature follows intensity'})
        saveThatShit(['cross corr combined bootstrapped' curvTypes{j}],p.OutputDirectory)
    end
end

%% ----------- Per- Branch Analysis -------- %%
    
%***********
%NOTE / WARNING - PER-BRANCH ANALYSIS IS NOT CURRENTLY NORMALIZED!
%....and hasn't been modified to suipport channel index input!!!!
%***********


allAvgInt = vertcat(avgBranchIntPerFrame{:});
allMaxInt = vertcat(maxBranchIntPerFrame{:});
allMedInt = vertcat(medBranchIntPerFrame{:});
allBranchLen = vertcat(branchLenPerFrame{:});
allBranchRad = vertcat(avgBranchRadiusPerFrame{:});        

for j = 1:nChan
    
    figure
    
    plot(allBranchRad,allAvgInt(:,j),'.','MarkerSize',25)
    hold on
    xlabel('Mean Branch Radius, nm')
    ylabel('Mean Intensity')
    title(['Intensity vs. Radius, Channel ' num2str(j)])
    
    figName = ['Combined branch radius vs. average intensity channel ' num2str(j)];
    print(pOpt{:},[outDir filesep figName '.eps']);
    print(pOptTif{:},[outDir filesep figName '.tif']);
    hgsave(gcf,[outDir filesep figName '.fig']);       
    
    figure
    
    plot(allBranchRad,allMaxInt(:,j),'.','MarkerSize',25)
    hold on
    xlabel('Mean Branch Radius, nm')
    ylabel('Max Intensity')
    title(['Intensity vs. Radius, Channel ' num2str(j)])
    
    figName = ['Combined branch radius vs. max intensity channel ' num2str(j)];
    print(pOpt{:},[outDir filesep figName '.eps']);
    print(pOptTif{:},[outDir filesep figName '.tif']);
    hgsave(gcf,[outDir filesep figName '.fig']);       
    
    figure
    
    plot(allBranchRad,allMedInt(:,j),'.','MarkerSize',25)
    hold on
    xlabel('Mean Branch Radius, nm')
    ylabel('Median Intensity')
    title(['Intensity vs. Radius, Channel ' num2str(j)])
    
    figName = ['Combined branch radius vs. median intensity channel ' num2str(j)];
    print(pOpt{:},[outDir filesep figName '.eps']);
    print(pOptTif{:},[outDir filesep figName '.tif']);
    hgsave(gcf,[outDir filesep figName '.fig']);       
    
end
   
%% ------------- Intensity Distribution / Normalization Analysis --------------- %%

movCols = lines(nMov);
minAll = Inf;
for k = iIntUse
    
    for l = 1:nChan
    
        currFig = fsFigure(.6);
        hold on

        for j = 1:nMov

            subplot(2,1,1)
            hold on
            [currN,currBin] = hist(allIntPerMov{j}(:,k,l),100);
            currN = currN ./ sum(currN);
            plot(currBin,currN,'color',movCols(j,:))
            xlabel([intNames{k} ', Channel ' num2str(p.ChannelIndex(l))])
            ylabel('Pixel Probability')
            legStr{j} = ['Movie ' num2str(j)];
            title(['Cortical ' intNames{k} ' sample distribution per-movie, Channel ' num2str(p.ChannelIndex(l))])    
            
            if p.NormalizeIntensities
                %Normalize to equal mean and standard deviation
                allIntPerMov{j}(:,k,l) = allIntPerMov{j}(:,k,l) - mean(allIntPerMov{j}(:,k,l));
                allIntPerMov{j}(:,k,l) = allIntPerMov{j}(:,k,l) ./ std(allIntPerMov{j}(:,k,l));
                minAll = min(min(allIntPerMov{j}(:,k,l)),minAll);
                
            end
            
        end
        for j = 1:nMov
                
            if p.NormalizeIntensities
                %Shift so that minimum of all movies is zero, because
                %people get weirded out if they see negative intensities,
                %even though they're fucking normalized
                allIntPerMov{j}(:,k,l) = allIntPerMov{j}(:,k,l) - minAll;
            
                subplot(2,1,2)
                hold on
                [currN,currBin] = hist(allIntPerMov{j}(:,k,l),100);
                currN = currN ./ sum(currN);
                plot(currBin,currN,'color',movCols(j,:))
                xlabel(['Normalized ' intNames{k} ', Channel ' num2str(p.ChannelIndex(l))])
                ylabel('Pixel Probability')
                title(['Normalized cortical ' intNames{k} ' sample distribution per-movie, Channel ' num2str(p.ChannelIndex(l))])    
                
            else
                subplot(2,1,1)
                hold on
                title('Normalization Disabled')                
                
            end
        end
        legend(legStr{:})        
        mfFigureExport(currFig,[p.OutputDirectory filesep 'Cortical ' intNames{k} ' sample distribution per-movie Channel ' num2str(p.ChannelIndex(l))])    
        
        
        
    end
end

%% ----------- Combined Intensity vs. Curvature Figures --------- %%
%Combines at the level of samples. Gives higher n, but will effectively
%weight cells by cortical area. Then again, if you're studying cortical
%structure... But cell-cell variability is higher??

nIntBins = 100;

allCurv = real(vertcat(allCurvPerMov{:}));
allInt = vertcat(allIntPerMov{:});

iGaussType = 1;
iMeanType = 2;
H = allCurv(:,iMeanType);
K = allCurv(:,iGaussType);

alpha = .05;
nBoot = 5e3;
bOpt.UseParallel = 'always';

curvMean = nan(nChan,nCurvTypes,numel(iIntUse),nIntBins);
curvSTD= nan(nChan,nCurvTypes,numel(iIntUse),nIntBins);
curvN = nan(nChan,nCurvTypes,numel(iIntUse),nIntBins);
curvSEM= nan(nChan,nCurvTypes,numel(iIntUse),nIntBins);
curvBootCI= nan(nChan,nCurvTypes,numel(iIntUse),nIntBins,2);

%Get the combined mean oversampling of the data
allNResamp = vertcat(nResampFacesPerFrame{:});
meanNResamp = mean(allNResamp);
curvSubSampCI= nan(nChan,nCurvTypes,numel(iIntUse),nIntBins,2);

for l = 1:nChan    
    
    for j = iIntUse 

        for k = 1:nCurvTypes
            
            % ---------- Average Curv Binned by Vs. Int ----- %%

            
            %intBins(j,:) = linspace(min(allInt(:,j,l)),max(allInt(:,j,l))+eps,nIntBins+1);
            pTiles(j,:) = linspace(0,100,nIntBins+1);
            intBins(j,:) = prctile(allInt(:,j),pTiles(j,:));
            ptCent(j,:) = pTiles(j,:) + [diff(pTiles(j,:)) / 2 0];%Percentiles centers for plotting
            intCent(j,:) = intBins(j,:) + [diff(intBins(j,:)) / 2 0];%Bin centers for plotting
            tic
            
            for m = 1:nIntBins
                currPts = allInt(:,j,l) >= intBins(j,m) & allInt(:,j,l) < intBins(j,m+1);
                currCurv = allCurv(currPts,k);
                curvMean(l,k,j,m) = mean(currCurv);
                curvSTD(l,k,j,m) = std(currCurv);
                curvN(l,k,j,m) = numel(currCurv);
                curvSEM(l,k,j,m) = curvSTD(l,k,j,m) / sqrt(curvN(l,k,j,m));
                
                %bsSamp = bootstrp(nBoot,@mean,currCurv,'Options' ,bOpt);
                if numel(currCurv) > 1
                    bsSamp = bootstrp(nBoot,@nanmean,currCurv);
                    curvBootCI(l,k,j,m,:) = prctile(bsSamp,[alpha/2 1-alpha/2]*100);       
                    ssSamp = subSampleBootStrap(nBoot,@nanmean,currCurv,round(numel(currCurv)/meanNResamp));
                    curvSubSampCI(l,k,j,m,:) = prctile(ssSamp,[alpha/2 1-alpha/2]*100);       
                end                                
                
            end
            toc
            
            % -------------- 2D Curve ----------- %
            
            currFig = figure;
            plot(squeeze(intCent(j,1:end-1)),squeeze(curvMean(l,k,j,:)))
            hold on
            plot(squeeze(intCent(j,1:end-1)),squeeze(curvBootCI(l,k,j,:,1)),'--')
            legend('Mean','95% C.I.')
            
            plot(squeeze(intCent(j,1:end-1)),squeeze(curvBootCI(l,k,j,:,2)),'--')
            xlim(intBins(j,[1 end-1]))
            plot(xlim,[0 0],'--r')
            xlabel([intNames{j} ' in channel ' num2str(l) ', a.u.'])
            ylabel([curvNames{k} ', ' curvUnits{k}])
            figName = [p.OutputDirectory filesep curvNames{k} ' versus ' intNames{j} ' channel ' num2str(l) ' plot'];
            mfFigureExport(currFig,figName)
            
            % -------------- 2D Curve by Percentile Regular BootStrap ----------- %
            
            currFig = figure;
            plot(squeeze(ptCent(j,1:end-1)),squeeze(curvMean(l,k,j,:)))
            hold on
            plot(squeeze(ptCent(j,1:end-1)),squeeze(curvBootCI(l,k,j,:,1)),'--')
            legend('Mean','95% C.I.')
            
            plot(squeeze(ptCent(j,1:end-1)),squeeze(curvBootCI(l,k,j,:,2)),'--')
            xlim(ptCent(j,[1 end-1]))
            plot(xlim,[0 0],'--r')
            xlabel(['Percentile of ' intNames{j} ' in channel ' num2str(l) ', a.u.'])
            ylabel([curvNames{k} ', ' curvUnits{k}])
            figName = [p.OutputDirectory filesep curvNames{k} ' versus ' intNames{j} ' channel ' num2str(l) ' plot by percentile'];
            mfFigureExport(currFig,figName)
            
            
             % -------------- 2D Curve by Percentile Sub-Sampled BootStrap ----------- %
            
            currFig = figure;
            plot(squeeze(ptCent(j,1:end-1)),squeeze(curvMean(l,k,j,:)))
            hold on
            plot(squeeze(ptCent(j,1:end-1)),squeeze(curvSubSampCI(l,k,j,:,1)),'--')
            legend('Mean','95% C.I.')
            
            plot(squeeze(ptCent(j,1:end-1)),squeeze(curvSubSampCI(l,k,j,:,2)),'--')
            xlim(ptCent(j,[1 end-1]))
            plot(xlim,[0 0],'--r')
            xlabel(['Percentile of ' intNames{j} ' in channel ' num2str(l) ', a.u.'])
            ylabel([curvNames{k} ', ' curvUnits{k}])
            figName = [p.OutputDirectory filesep curvNames{k} ' versus ' intNames{j} ' channel ' num2str(l) ' plot by percentile subsampled bootstrap'];
            
                        
            % -------------- 2D Histogram ----------- %
            
            currFig = figure;                        
            [N,C] = hist3([allInt(:,j,l),allCurv(:,k)],[200 200]);
            imagesc(C{1},C{2},log10(N')),axis xy
            hold on
            plot(squeeze(intCent(j,1:end-1)),squeeze(curvMean(l,k,j,:)))
            hold on
            plot(squeeze(intCent(j,1:end-1)),squeeze(curvBootCI(l,k,j,:,1)),'--')
            legend('Mean','95% C.I.')
            plot(squeeze(intCent(j,1:end-1)),squeeze(curvBootCI(l,k,j,:,2)),'--')
            xlabel([intNames{j} ' in channel ' num2str(l) ', a.u.'])
            ylabel([curvNames{k} ', ' curvUnits{k}])
            plot(xlim,[0 0],'--r')
            colorbar
            title({'Curvature vs. intensity 2D Histogram',...
                    'Color indicates log10 of sample count',...
                    ['n=' num2str(nMov) ' cells, ' num2str(sum(nFramesPerMov)) ' time points, ' num2str(size(allInt,1)) ' samples']})                            
            useRange = [ 0 log10(prctile(N(isfinite(N(:)) & N(:) > 0),99))];
            
            caxis(useRange)
            colormap gray
            figName = [p.OutputDirectory filesep curvNames{k} ' versus ' intNames{j} ' channel ' num2str(l) ' 2D Histogram'];
            mfFigureExport(currFig,figName)
            
        end
        
        % ------------- Curv Category Vs. Intensity ------- %
        

        nPerCat = zeros(nCurvCat,nIntBins-1);
        
        
        for m = 1:nIntBins                        
            currPts = allInt(:,j,l) >= intBins(j,m) & allInt(:,j,l) < intBins(j,m+1);            
            for iCat = 1:nCurvCat
                currPtsCat = (K > curvCatKLims(iCat,1) & K <= curvCatKLims(iCat,2) & H > curvCatHLims(iCat,1) & H <= curvCatHLims(iCat,2)) & currPts;                
                nPerCat(iCat,m) = nnz(currPtsCat);                                        
            end            
        end

            
        % --- Total Sample Count Per Category --- %
    
        currFig = figure;    
        %bHan = bar(intCent(j,1:end-1),nPerCat',1,'stacked');
        bHan = area(intCent(j,1:end-1),nPerCat');
        xlabel([intNames{j} ' in channel ' num2str(l) ', a.u.'])
        ylabel('Number of Surface Samples')
        legend(curvCatNames,'Location','NorthEastOutside')
        for iCat = 1:nCurvCat
            set(bHan(iCat),'FaceColor',curvCatColors(iCat,:))
            set(bHan(iCat),'EdgeColor',curvCatColors(iCat,:))
        end
        xlim([min(intCent(j,:)),max(intCent(j,1:end-1))])
        figName = [p.OutputDirectory filesep 'Curvature Category Versus ' intNames{j} ' Channel ' num2str(l) ' sample count'];
        mfFigureExport(currFig,figName)

        % --- Fraction of Samples Per Category --- %

        nTotInCat = sum(nPerCat,1);

        fracPerCat = bsxfun(@rdivide,nPerCat,nTotInCat);
        %fracPerCat = bsxfun(@divide nPerCat 

        currFig = figure;    
         %bHan = bar(intCent(j,1:end-1),fracPerCat',1,'stacked');
        bHan = area(intCent(j,1:end-1),fracPerCat');
        xlabel([intNames{j} ' in channel ' num2str(l) ', a.u.'])
        ylabel('Fraction of Surface Samples')
        legend(curvCatNames,'Location','NorthEastOutside')
        for iCat = 1:nCurvCat
            set(bHan(iCat),'FaceColor',curvCatColors(iCat,:))
            set(bHan(iCat),'EdgeColor',curvCatColors(iCat,:))
        end
        xlim([min(intCent(j,:)),max(intCent(j,1:end-1))])
        figName = [p.OutputDirectory filesep 'Curvature Category Versus ' intNames{j} ' Channel ' num2str(l) ' sample fraction'];
        mfFigureExport(currFig,figName)    
                    
        % --- Fraction of Samples Per Category by Percentile --- %

        nTotInCat = sum(nPerCat,1);

        fracPerCat = bsxfun(@rdivide,nPerCat,nTotInCat);
        %fracPerCat = bsxfun(@divide nPerCat 

        currFig = figure;    
         %bHan = bar(intCent(j,1:end-1),fracPerCat',1,'stacked');
        bHan = area(ptCent(j,1:end-1),fracPerCat');
        xlabel(['Percentile of ' intNames{j} ' in channel ' num2str(l) ', a.u.'])
        ylabel('Fraction of Surface Samples')
        legend(curvCatNames,'Location','NorthEastOutside')
        for iCat = 1:nCurvCat
            set(bHan(iCat),'FaceColor',curvCatColors(iCat,:))
            set(bHan(iCat),'EdgeColor',curvCatColors(iCat,:))
        end
        xlim([min(ptCent(j,:)),max(ptCent(j,1:end-1))])
        figName = [p.OutputDirectory filesep 'Curvature Category Versus ' intNames{j} ' Channel ' num2str(l) ' sample fraction percentile'];
        mfFigureExport(currFig,figName)    

        
    end
end

%% ----------- Per-Cell Curvature and Intensity Correlation Overlays ----------- %%  

cellCols = rand(nMov,3);

for l = 1:nChan
    
    for k = 1:nCurvTypes

        for j = iIntUse
            
            % ---------- Average Curv Binned by Vs. Int ----- %%
            
            currFig = figure;
            hold on
            for n = 1:nMov
                
                for m = 1:nIntBins
                    currPts = allIntPerMov{n}(:,j,l) >= intBins(j,m) & allIntPerMov{n}(:,j,l) < intBins(j,m+1);
                    currCurv = real(allCurvPerMov{n}(currPts,k));
                    curvMeanPer(n,l,k,j,m) = mean(currCurv);
                    curvSTDPer(n,l,k,j,m) = std(currCurv);
                    curvNPer(n,l,k,j,m) = numel(currCurv);
                    curvSEMPer(n,l,k,j,m) = curvSTDPer(n,l,k,j,m) / sqrt(curvNPer(n,l,k,j,m));
                end

                plot(squeeze(intCent(j,1:end-1)),squeeze(curvMeanPer(n,l,k,j,:)),'color',cellCols(n,:))
%Gets too confusing with uncertainties displayed                
                 %                 hold on
%                 plot(squeeze(intCent(j,1:end-1)),squeeze(curvMeanPer(n,l,k,j,:)) + 1.96*squeeze(curvSEMPer(n,l,k,j,:)),'--','color',cellCols(n,:))
                
%                plot(squeeze(intCent(j,1:end-1)),squeeze(curvMeanPer(n,l,k,j,:)) - 1.96*squeeze(curvSEMPer(n,l,k,j,:)),'--','color',cellCols(n,:))
                
                
            end
            legend(arrayfun(@num2str,1:nMov,'Unif',0))
            xlim(intBins(j,[1 end-1]))
            plot(xlim,[0 0],'--r')
            xlabel([intNames{j} ' in channel ' num2str(l) ', a.u.'])
            ylabel([curvNames{k} ', ' curvUnits{k}])
            figName = [p.OutputDirectory filesep curvNames{k} ' versus ' intNames{j} ' channel ' num2str(l) ' per-cell curv overlay plot'];
            title('Overlay of per-cell curvature correlation curves')
            mfFigureExport(currFig,figName)
            
        end

    end
end

%% ----------- Combined Average of Per-Cell Correlations ------- %%

%Gets combined correlations treating each cell as one datapoint. This puts
%greater emphasis on cell-cell variability, weights equally per cell, but
%can therefore give less weight per cortical area in larger cells.

%%NOTE PER CELL IS NOT BOOTSTRAPPED YET

if nMov > 1
    
    curvMeanOfPerCell = nan(nChan,nCurvTypes,numel(iIntUse),nIntBins);
    curvSTDOfPerCell = nan(nChan,nCurvTypes,numel(iIntUse),nIntBins);
    curvNOfPerCell = nan(nChan,nCurvTypes,numel(iIntUse),nIntBins);
    curvSEMOfPerCell = nan(nChan,nCurvTypes,numel(iIntUse),nIntBins);
    curvBootCIOfPerCell = nan(nChan,nCurvTypes,numel(iIntUse),nIntBins,2);

    for l = 1:nChan    

        for k = 1:nCurvTypes

            for j = iIntUse

                curvMeanOfPerCell(l,k,j,:) = squeeze(nanmean(curvMeanPer(:,l,k,j,:),1));
                curvSTDOfPerCell(l,k,j,:) = squeeze(nanstd(curvMeanPer(:,l,k,j,:),[],1));
                curvNOfPerCell(l,k,j,:) = squeeze(sum(~isnan(curvMeanPer(:,l,k,j,:)),1));
                curvSEMOfPerCell(l,k,j,:) = curvSTDOfPerCell(l,k,j,:) ./ sqrt(curvNOfPerCell(l,k,j,:));

    %             
                for m = 1:nIntBins
                    bsSamp = bootstrp(nBoot,@nanmean,squeeze(curvMeanPer(:,l,k,j,m)));
                    curvBootCIOfPerCell(l,k,j,m,:) = prctile(bsSamp,[alpha/2 1-alpha/2]*100);       
                end     

    %             curvBootCIOfPerCell(l,k,j,:) = 
    %             
                % ---------- Average of Per-Cell Bootstrapped ----- %

                currFig = figure;

                plot(intCent(j,1:end-1),squeeze(curvMeanOfPerCell(l,k,j,:))')            
                hold on
                plot(intCent(j,1:end-1),squeeze(curvBootCIOfPerCell(l,k,j,:,1))','--')
                legend('Mean','95% C.I.')
                plot(intCent(j,1:end-1),squeeze(curvBootCIOfPerCell(l,k,j,:,2))','--')
                xlim(intBins(j,[1 end-1]))
                plot(xlim,[0 0],'--r')
                xlabel([intNames{j} ' in channel ' num2str(l) ', a.u.'])
                ylabel([curvNames{k} ', ' curvUnits{k}])
                figName = [p.OutputDirectory filesep curvNames{k} ' versus ' intNames{j} ' channel ' num2str(l) ' average of per-cell curv correlations plot'];
                title({'Average of per-cell correlations',...
                    ['n=' num2str(nMov) ' cells, ' num2str(sum(nFramesPerMov)) ' time points']})
                mfFigureExport(currFig,figName)

                % ---------- Average of Per-Cell Bootstrapped and by percentile ----- %

                currFig = figure;

                plot(ptCent(j,1:end-1),squeeze(curvMeanOfPerCell(l,k,j,:))')            
                hold on
                plot(ptCent(j,1:end-1),squeeze(curvBootCIOfPerCell(l,k,j,:,1))','--')
                legend('Mean','95% C.I.')
                plot(ptCent(j,1:end-1),squeeze(curvBootCIOfPerCell(l,k,j,:,2))','--')
                xlim(ptCent(j,[1 end-1]))
                plot(xlim,[0 0],'--r')
                xlabel(['Percentile of ' intNames{j} ' in channel ' num2str(l) ', a.u.'])
                ylabel([curvNames{k} ', ' curvUnits{k}])
                figName = [p.OutputDirectory filesep curvNames{k} ' versus ' intNames{j} ' channel ' num2str(l) ' average of per-cell curv correlations plot by percentile'];
                title({'Average of per-cell correlations',...
                    ['n=' num2str(nMov) ' cells, ' num2str(sum(nFramesPerMov)) ' time points']})
                mfFigureExport(currFig,figName)

            end

        end

    end
else
    curvMeanOfPerCell = nan;
    curvSTDOfPerCell = nan;
    curvNOfPerCell = nan;
    curvSEMOfPerCell = nan;
    curvBootCIOfPerCell = nan;
end

%% ------------------ Output ------------ %%

outVars = {'curvMean','curvSTD','curvN','curvSEM','intTypes','intNames','iIntUse','chanPaths','curvMeanOfPerCell','curvSTDOfPerCell','curvNOfPerCell','curvSEMOfPerCell',...
           'nPerCat'};

save([p.OutputDirectory filesep 'combined analysis.mat'],outVars{:})
