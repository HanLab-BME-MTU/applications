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
nFramesPer = nan(nMov,1);
nChanPer = nan(nMov,1);
pixelSizePer = nan(nMov,1);

curvConv = cell(nMov,1);
[intTypes intNames] = getIntTypeFields;
nIntTypes = numel(intTypes);


%% ------------- Per-Movie Loading and Processing --------------- %%


for iMov = 1:nMov
    
    nFramesPer(iMov) = MA(iMov).nFrames_;
    pixelSizePer(iMov) = MA(iMov).pixelSize_;
    nChanPer(iMov) = numel(MA(iMov).channels_);
    
    iMAProc = MA(iMov).getProcessIndex('MaskedIntensity3DProcess',1 ,0);        
    assert(~isempty(iMAProc),['Movie ' num2str(iMov) ' has no intensity analysis!']);    
        
    intAn{iMov} = MA(iMov).processes_{iMAProc}.loadChannelOutput(iProcChan);
    
    [curvTypes,curvNames,curvUnits,curvConv{iMov}] = getCurveTypeFields(pixelSizePer(iMov),useMicrons);
    nCurvTypes = numel(curvTypes);               
        
    if nFramesPer(iMov) > 1
    
        %TEEEEMPPP ADD MULTI-CHANNEL SUPPORT??? Probably wont see
        %multi-channel timelaps for years....
        ccPerMov(iMov,:,:) = intAn{iMov}.ccIntCurv(1,:,:);
        tIntPer(iMov) = MA(iMov).timeInterval_;                
    else
        tIntPer(iMov) = 0;
    end
        
    
    %Extract and combine per-frame data
                
    
    nSurfPtsPerFrame{iMov} = cellfun(@numel,{intAn{iMov}.branchProfiles(:).surfPixInd});
    nSurfPotsTot(iMov) = sum(nSurfPtsPerFrame{iMov});
    allCurvPerMov{iMov} = nan(nSurfPotsTot(iMov),nCurvTypes);    
    for j = 1:nCurvTypes
        allCurvPerMov{iMov}(:,j) = vertcat( intAn{iMov}.branchProfiles(:).(curvTypes{j}) ) .* curvConv{iMov}(j);            
    end
    
    allIntPerMov{iMov} = nan(nSurfPotsTot(iMov),nIntTypes,nChanPer(iMov));    
    for j = 1:nIntTypes    
        allIntPerMov{iMov}(:,j,:) = vertcat(intAn{iMov}.branchProfiles(:).(intTypes{j}));     
    end
    
    for iFrame = 1:nFramesPer(iMov)
        
        nBranchPerFrame{iMov}(iFrame) = numel(intAn{iMov}.branchProfiles(iFrame).branchTipPixelLenVsDepth);            
        hasLvD = cellfun(@(x)(~isempty(x) && ~any(isnan(x(:)))),intAn{iMov}.branchProfiles(iFrame).branchTipPixelLenVsDepth);                        
        
        avgBranchRadiusPerFrame{iMov,iFrame} = cellfun(@(x)(mean(x(:,2))),intAn{iMov}.branchProfiles(iFrame).branchTipPixelLenVsDepth(hasLvD)) .* MA(iMov).pixelSize_;        
        
        for k = 1:nChanPer(iMov);
            
            avgBranchIntPerFrame{iMov,iFrame}(:,k) = cellfun(@(x)(mean(x(:,k))),intAn{iMov}.branchProfiles(iFrame).branchTipPixelInt(hasLvD));
            maxBranchIntPerFrame{iMov,iFrame}(:,k) = cellfun(@(x)(max(x(:,k))),intAn{iMov}.branchProfiles(iFrame).branchTipPixelInt(hasLvD));
            medBranchIntPerFrame{iMov,iFrame}(:,k) = cellfun(@(x)(median(x(:,k))),intAn{iMov}.branchProfiles(iFrame).branchTipPixelInt(hasLvD));

        end                
        branchLenPerFrame{iMov,iFrame} = cellfun(@(x)(max(x(:,1))),intAn{iMov}.branchProfiles(iFrame).branchTipPixelLenVsDepth(hasLvD)) .* MA(iMov).pixelSize_;
    end
end

assert(numel(unique(nChanPer))==1,'All movies must have the same number of channels!!')
nChan = nChanPer(1);


if numel(unique(tIntPer)) > 1 || numel(unique(nFramesPer)) > 1
    %Obviously this isn't strictly true, we're just using this constraint
    %for now.
    warning('All movies should have the same time-interval and number of frames to run combined cross-corr!!')    
end


nCurvType = numel(curvTypes);


%% ------------- Cross-Correlation Analysis  --------------- %%

tInt = tIntPer(1);    
nFrames = nFramesPer(1);
isTLapse = tInt ~= 0;

if isTLapse

    maxDelay = (size(ccPerMov,3)-1)/2;
    tLags = -tInt*maxDelay:tInt:tInt*maxDelay;

    for j = 1:nCurvType

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
   
%% ------------- Combined Curvature Analysis --------------- %%


jk=1;