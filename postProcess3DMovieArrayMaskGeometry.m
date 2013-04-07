function postProcess3DMovieArrayMaskGeometry(MA,varargin)
%POSTPROCESS3DMOVIEARRAYMASKGEOMETRY calculates various statistics regarding the geometry of the masks in the input movie array 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% UNDER CONSTRUCTION
% 
% 
% 
% 

%% ------------- Parameters -------------- %%

binSz = 1/1e4;%Bin size for curvature histograms, in 1/nm
gcBinSz = binSz^2*10; %Gaussian curvature is product of principle curvatures
nPixMaxCurv = 2;%Radius in pixels of the maximum-curvature-radius to set histogram limits at

perMovDirName = 'mask geometry post processing';%Directory for saving individual move results in movie output directory

%Print options for saving figures to disk
pOpt = {'-r300',...% dpi = 300
        '-depsc2'};% use eps format

%% ------------- Input ------------ %%


%Parse all the inputs
ip = inputParser;
ip.FunctionName = mfilename;
ip.addRequired('MA',@(x)(isa(x,'MovieData3D')));
ip.addParamValue('ChannelIndex',1,@(x)(numel(x) == 1 && isposint(x)));
ip.addParamValue('SampRad',2e3,@(x)(numel(x) == 1 && x > 0));
ip.addParamValue('BatchMode',false,(@(x)(numel(x)==1)));
ip.addParamValue('OutputDirectory',[],(@(x)(exist(x,'dir'))));
ip.parse(MA,varargin{:});

p = ip.Results;

if isempty(p.OutputDirectory)
    p.OutputDirectory = uigetdir(pwd,'Select a directory to store the results:');
    if p.OutputDirectory == 0
        error('You must select an output directory!')
    end
end

%% ------------ Init --------------- %%


nMovies = numel(MA);
iMGProc = zeros(nMovies,1);
hasMG = false(nMovies,1);
nFramesPerMov = nan(nMovies,1);
pixSizePerMov = nan(nMovies,1);
histBinsPerMov = cell(nMovies,1);
gcHistBinsPerMov = cell(nMovies,1);
dpcHistBinsPerMov = cell(nMovies,1);
mapcHistBinsPerMov = cell(nMovies,1);
gcHistPerMov = cell(nMovies,1);%Gauss curv hist
mcHistPerMov = cell(nMovies,1);%mean curv hist
dpcHistPerMov = cell(nMovies,1);%Diff of principle curv hist
mapcHistPerMov = cell(nMovies,1);%Maximum absolute principle curv hist
nBinsPerMov = nan(nMovies,1);
ngcBinsPerMov = nan(nMovies,1);
ndpcBinsPerMov = nan(nMovies,1);
nmapcBinsPerMov = nan(nMovies,1);

physUnits = false(nMovies,1);
%timIntPerMov = nan(nMovies,1);

gcFigPM = zeros(nMovies,1);%gauss curv fig handles per movie
mcFigPM = zeros(nMovies,1);%mean curv fig handles per movie
dpcFigPM = zeros(nMovies,1);%diff of princip curv fig handles per movie

if ~p.BatchMode
    wtBar = waitbar(0,'Please wait, performing post-processing on each movie...');
end        


%% -------------- Per-Movie Processing ----------------- %%





for iMov = 1:nMovies
    
    tmp = MA(iMov).getProcessIndex('MaskGeometry3DProcess',1,~p.BatchMode);
    
    if ~isempty(tmp) && MA(iMov).processes_{tmp}.checkChannelOutput(p.ChannelIndex);
                        
        
        iMGProc(iMov) = tmp;
        hasMG(iMov) = true;        
        if isempty(MA(iMov).eventTimes_)
            nFramesPerMov(iMov) = MA(iMov).nFrames_;
        else
            %TEMP - this doesn't currently support a start frame > 1!!!
            nFramesPerMov(iMov) = MA(iMov).eventTimes_(2);
        end        
        pixSizePerMov(iMov) = MA(iMov).pixelSize_;
        %timeIntPerMov(iMov) = MA(iMov).timeInterval_;
        
        %Set up the curvature histogram bins for this movie based on the
        %pixel size - no object can have higher curvature than a single
        %pixel.
        %%% TEMP !!!!!!!!!!!!!!!! 
%         histBinsPerMov{iMov} = -1/(pixSizePerMov(iMov)*nPixMaxCurv):binSz:1/(pixSizePerMov(iMov)*nPixMaxCurv);
%         nBinsPerMov(iMov) = numel(histBinsPerMov{iMov});
%         gcHistBinsPerMov{iMov} = -(1/pixSizePerMov(iMov)/nPixMaxCurv)^2:gcBinSz:(1/pixSizePerMov(iMov)/nPixMaxCurv)^2;
%         ngcBinsPerMov(iMov) = numel(gcHistBinsPerMov{iMov});
%         dpcHistBinsPerMov{iMov} = 0:binSz:1/(pixSizePerMov(iMov)*nPixMaxCurv);
%         ndpcBinsPerMov(iMov) = numel(dpcHistBinsPerMov{iMov});
        %Use hard-coded bins so we can easily
        %compare across movies/conditions etc..... i know, i know...
        histBinsPerMov{iMov} = -2e-3:binSz:2e-3;
        nBinsPerMov(iMov) = numel(histBinsPerMov{iMov});
        gcHistBinsPerMov{iMov} = -2e-6:gcBinSz:2e-6;
        ngcBinsPerMov(iMov) = numel(gcHistBinsPerMov{iMov});
        dpcHistBinsPerMov{iMov} = 0:binSz:2e-3;
        ndpcBinsPerMov(iMov) = numel(dpcHistBinsPerMov{iMov});
        mapcHistBinsPerMov{iMov} = 0:binSz:5e-3;
        nmapcBinsPerMov(iMov) = numel(mapcHistBinsPerMov{iMov});
        
        gcHistPerMov{iMov} = nan(nFramesPerMov(iMov),ngcBinsPerMov(iMov));
        mcHistPerMov{iMov} = nan(nFramesPerMov(iMov),nBinsPerMov(iMov));           
        dpcHistPerMov{iMov} = nan(nFramesPerMov(iMov),ndpcBinsPerMov(iMov));
        mapcHistPerMov{iMov} = nan(nFramesPerMov(iMov),nmapcBinsPerMov(iMov));       
        
        currOutDir = [MA(iMov).outputDirectory_ filesep perMovDirName];
        mkClrDir(currOutDir);
        
        physUnits(iMov) = MA(iMov).processes_{iMGProc(iMov)}.funParams_.PhysicalUnits;
        
        
        for iFrame = 1:nFramesPerMov(iMov)
        
            currMaskGeo = MA(iMov).processes_{iMGProc(iMov)}.loadChannelOutput(p.ChannelIndex,iFrame);
            
            if ~physUnits(iMov)
                %Scale these to physical units of 1/nm. We only scale by
                %the XY pixel size because the voxels have been made
                %symmetric
                currMaskGeo.GaussianCurvature = currMaskGeo.GaussianCurvature ./ (pixSizePerMov(iMov)^2);
                currMaskGeo.MeanCurvature = currMaskGeo.MeanCurvature ./ pixSizePerMov(iMov);
                currMaskGeo.CurvaturePC1 = currMaskGeo.CurvaturePC1 ./ pixSizePerMov(iMov);
                currMaskGeo.CurvaturePC2 = currMaskGeo.CurvaturePC2 ./ pixSizePerMov(iMov);
            end
            
            %Perform local averaging on curvature data to reduce noise and
            %match with image sampling data. 
            %Even though the curv data has been converted to physical units
            %the 3d mesh is still in pixel coord so we convert the sampling
            %radius here.
            locAvgMaskGeo = calcLocalAvgCurvatures(currMaskGeo,p.SampRad / pixSizePerMov(iMov),1,1);
            
            
            %Get normalized curvature histograms
            gcHistPerMov{iMov}(iFrame,:) = histc(locAvgMaskGeo.LocMeanGaussianCurvature,gcHistBinsPerMov{iMov});
            gcHistPerMov{iMov}(iFrame,:) = gcHistPerMov{iMov}(iFrame,:) ./ sum(gcHistPerMov{iMov}(iFrame,:));
            mcHistPerMov{iMov}(iFrame,:) = histc(locAvgMaskGeo.LocMeanMeanCurvature,histBinsPerMov{iMov});
            mcHistPerMov{iMov}(iFrame,:) = mcHistPerMov{iMov}(iFrame,:) ./ sum(mcHistPerMov{iMov}(iFrame,:));
            dpcHistPerMov{iMov}(iFrame,:) = histc(real(locAvgMaskGeo.LocMeanCurvaturePC1) - real(locAvgMaskGeo.LocMeanCurvaturePC2),dpcHistBinsPerMov{iMov});
            dpcHistPerMov{iMov}(iFrame,:) = dpcHistPerMov{iMov}(iFrame,:) ./ sum(dpcHistPerMov{iMov}(iFrame,:));                    
            mapcHistPerMov{iMov}(iFrame,:) = histc(locAvgMaskGeo.LocMeanMaxAbsCurvature,mapcHistBinsPerMov{iMov});
            mapcHistPerMov{iMov}(iFrame,:) = mapcHistPerMov{iMov}(iFrame,:) ./ sum(mapcHistPerMov{iMov}(iFrame,:));
            
            %Get per-frame simple stats. TEMP - preallocate arrays!
            gcMeanPerMov{iMov}(iFrame) = mean(locAvgMaskGeo.LocMeanGaussianCurvature);
            gcMedPerMov{iMov}(iFrame) = median(locAvgMaskGeo.LocMeanGaussianCurvature);
            gcVarPerMov{iMov}(iFrame) = var(locAvgMaskGeo.LocMeanGaussianCurvature);
            mcMeanPerMov{iMov}(iFrame) = mean(locAvgMaskGeo.LocMeanMeanCurvature);
            mcMedPerMov{iMov}(iFrame) = median(locAvgMaskGeo.LocMeanMeanCurvature);
            mcVarPerMov{iMov}(iFrame) = var(locAvgMaskGeo.LocMeanMeanCurvature);
            mapcMeanPerMov{iMov}(iFrame) = mean(locAvgMaskGeo.LocMeanMaxAbsCurvature);
            mapcMedPerMov{iMov}(iFrame) = median(locAvgMaskGeo.LocMeanMaxAbsCurvature);
            mapcVarPerMov{iMov}(iFrame) = var(locAvgMaskGeo.LocMeanMaxAbsCurvature);
            [~,iMode] = max(mapcHistPerMov{iMov}(iFrame,:));%This is a sloppy as hell way to do this, but quick and dirty will have to do.
            mapcHistModePerMov{iMov}(iFrame) = mapcHistBinsPerMov{iMov}(iMode);
            
            %Store number of curv points
            nCurvPtsPerMov{iMov}(iFrame) = numel(locAvgMaskGeo.LocMeanGaussianCurvature);            
            
        end
                
        
        % ---- Gauss Curv Figure ---- %%
        
        if p.BatchMode
            gcFigPM(iMov) = figure('Visible','off');
        else
            gcFigPM(iMov) = figure;
        end
        set(gcFigPM(iMov),'DefaultAxesColorOrder',jet(nFramesPerMov(iMov)));
        %frameTimes = 0:timeIntPerMov(iMov):(timeIntPerMov(iMov)*(nFramesPerMov(iMov)-1));
        %imagesc(gcHistBinsPerMov{iMov},frameTimes,gcHistPerMov{iMov})
        plot(repmat(gcHistBinsPerMov{iMov},nFramesPerMov(iMov),1)',gcHistPerMov{iMov}')
        xlabel('Gaussian Curvature, 1/nm^2')
        ylabel('Normalized Histogram')
        title('Gaussian Curvature Distribution Over Time (blue=first frame,red=last frame)')
        print(pOpt{:},[currOutDir filesep 'Gaussian Curvature Distribution.eps']);
        
        
        % ---- Mean Curv Figure ---- %%
        
        if p.BatchMode
            mcFigPM(iMov) = figure('Visible','off');
        else
            mcFigPM(iMov) = figure;
        end
        set(mcFigPM(iMov),'DefaultAxesColorOrder',jet(nFramesPerMov(iMov)));
        %frameTimes = 0:timeIntPerMov(iMov):(timeIntPerMov(iMov)*(nFramesPerMov(iMov)-1));
        %imagesc(gcHistBinsPerMov{iMov},frameTimes,gcHistPerMov{iMov})
        plot(repmat(histBinsPerMov{iMov},nFramesPerMov(iMov),1)',mcHistPerMov{iMov}')
        xlabel('Mean Curvature, 1/nm')
        ylabel('Normalized Histogram')
        title('Mean Curvature Distribution Over Time (blue=first frame,red=last frame)')
        print(pOpt{:},[currOutDir filesep 'Gaussian Curvature Distribution.eps']);
        
        
        % ---- Diff of PC Figure ---- %%
        
        if p.BatchMode
            dpcFigPM(iMov) = figure('Visible','off');
        else
            dpcFigPM(iMov) = figure;
        end
        set(dpcFigPM(iMov),'DefaultAxesColorOrder',jet(nFramesPerMov(iMov)));
        %frameTimes = 0:timeIntPerMov(iMov):(timeIntPerMov(iMov)*(nFramesPerMov(iMov)-1));
        %imagesc(gcHistBinsPerMov{iMov},frameTimes,gcHistPerMov{iMov})
        plot(repmat(dpcHistBinsPerMov{iMov},nFramesPerMov(iMov),1)',dpcHistPerMov{iMov}')
        xlabel('Difference of Principle Curvatures, 1/nm')
        ylabel('Normalized Histogram')
        title('Diff. of Principle Curvatures Distribution Over Time (blue=first frame,red=last frame)')
        print(pOpt{:},[currOutDir filesep 'Difference of Principle Curvature Distribution.eps']);
        
         % ---- Max Abs PC Figure ---- %%
        
        if p.BatchMode
            mapcFigPM(iMov) = figure('Visible','off');
        else
            mapcFigPM(iMov) = figure;
        end
        set(mapcFigPM(iMov),'DefaultAxesColorOrder',jet(nFramesPerMov(iMov)));
        %frameTimes = 0:timeIntPerMov(iMov):(timeIntPerMov(iMov)*(nFramesPerMov(iMov)-1));
        %imagesc(gcHistBinsPerMov{iMov},frameTimes,gcHistPerMov{iMov})
        plot(repmat(mapcHistBinsPerMov{iMov},nFramesPerMov(iMov),1)',mapcHistPerMov{iMov}')
        xlabel('Max Absolute Principle Curvature, 1/nm')
        ylabel('Normalized Histogram')
        title('Max Absolute Principle Curvature Distribution Over Time (blue=first frame,red=last frame)')
        print(pOpt{:},[currOutDir filesep 'Max Absolute Principle Curvature Distribution.eps']);
        
        if p.BatchMode
            close(mcFigPM(iMov));
            close(gcFigPM(iMov));
        end
                
        
    else
        warning('MIGRATION3D:ppMaskGeo:noMaskGeo',...
            ['Movie ' num2str(iMov) ' does not have valid mask geometry analysis - not analyzing!']);
    end
    
    if ~p.BatchMode        
        waitbar(iMov/nMovies,wtBar)
    end
    
    
end
    
  
    
if ~p.BatchMode && ishandle(wtBar)    
    close(wtBar);
end



%% ------------ Combined Analysis of All Movies -------------- %%

compFig = figure;
hold on
movColors = jet(nMovies);

for iMov = 1:nMovies
    
    figure

    meanDiffHist(iMov,:) = nanmean(dpcHistPerMov{iMov},1);
    stdDiffHist(iMov,:) = nanstd(dpcHistPerMov{iMov},[],1);    
    
    plotTransparent(dpcHistBinsPerMov{iMov},meanDiffHist(iMov,:),stdDiffHist(iMov,:),movColors(iMov,:),.5,0);        
    hold on
    plot(dpcHistBinsPerMov{iMov},meanDiffHist(iMov,:),'color',movColors(iMov,:))
    xlim([1e-4 1.9e-3])
    xlabel('Difference of Principle Curvatures, 1/nm')
    ylabel('Probability')
    title('Diff. of Principle Curvatures Distribution Averaged Over All Frames')
    currOutDir = [MA(iMov).outputDirectory_ filesep perMovDirName];
    saveThatShit([currOutDir filesep 'Averaged Difference of Principle Curvature Distribution'],currOutDir);
    
    figure
    
    meanMeanHist(iMov,:) = nanmean(mcHistPerMov{iMov},1);
    stdMeanHist(iMov,:) = nanstd(mcHistPerMov{iMov},[],1);    
    
    plotTransparent(histBinsPerMov{iMov},meanMeanHist(iMov,:),stdMeanHist(iMov,:),movColors(iMov,:),.5,0);        
    hold on
    plot(histBinsPerMov{iMov},meanMeanHist(iMov,:),'color',movColors(iMov,:))
    xlim([-2e-3 2e-3])
    xlabel('Mean Curvature, 1/nm')
    ylabel('Probability')
    title('Mean Curvature Distribution Averaged Over All Frames')
    currOutDir = [MA(iMov).outputDirectory_ filesep perMovDirName];
    saveThatShit([currOutDir filesep 'Averaged Mean Curvature Distribution'],currOutDir);
    
    figure
    
    meanGaussHist(iMov,:) = nanmean(gcHistPerMov{iMov},1);
    stdGaussHist(iMov,:) = nanstd(gcHistPerMov{iMov},[],1);    
    
    plotTransparent(gcHistBinsPerMov{iMov},meanGaussHist(iMov,:),stdGaussHist(iMov,:),movColors(iMov,:),.5,0);        
    hold on
    plot(gcHistBinsPerMov{iMov},meanGaussHist(iMov,:),'color',movColors(iMov,:));        
    xlim([-2e-6 2e-6])
    xlabel('Gaussian Curvature, 1/nm^2')
    ylabel('Probability')
    title('Gaussian Curvature Distribution Averaged Over All Frames')
    currOutDir = [MA(iMov).outputDirectory_ filesep perMovDirName];
    saveThatShit([currOutDir filesep 'Averaged Gaussian Curvature Distribution'],currOutDir);
    
    
    figure
    
    meanMAPCHist(iMov,:) = nanmean(mapcHistPerMov{iMov},1);
    stdMAPCHist(iMov,:) = nanstd(mapcHistPerMov{iMov},[],1);    
    
    plotTransparent(mapcHistBinsPerMov{iMov},meanMAPCHist(iMov,:),stdMAPCHist(iMov,:),movColors(iMov,:),.5,0);        
    hold on
    plot(mapcHistBinsPerMov{iMov},meanMAPCHist(iMov,:),'color',movColors(iMov,:));        
    xlim([0 5e-3])
    xlabel('Max. Abs. Curvature Cmponent, 1/nm')
    ylabel('Probability')
    title('Maximum Absolute Curvature Component Distribution Averaged Over All Frames')
    currOutDir = [MA(iMov).outputDirectory_ filesep perMovDirName];
    saveThatShit([currOutDir filesep 'Averaged Max Abs Curv Component Distribution'],currOutDir);    
    
    
end

%--------- All Movie Distribution Averages -------- %

figure

combOutDir = p.OutputDirectory;
combMeanHist = mean(meanMeanHist,1);
combMeanHistSTD = std(meanMeanHist,[],1);
%TEEEEMMP TEMP TEMP !!! FUCKING THESIS !!!
plotTransparent(histBinsPerMov{1},combMeanHist,combMeanHistSTD,[0 0 1],.5,0);        
hold on
plot(histBinsPerMov{1},combMeanHist)
xlim([-2e-3 2e-3])
xlabel('Mean Curvature, 1/nm')
ylabel('Probability')
title(['Mean Curvature Distribution Averaged Over All Movies, n =' num2str(nMovies)])
saveThatShit(['Combined Mean Curvature Distribution'],combOutDir);

figure

combGaussHist = mean(meanGaussHist,1);
combGaussHistSTD = std(meanGaussHist,[],1);

plotTransparent(gcHistBinsPerMov{1},combGaussHist,combGaussHistSTD,[0 0 1],.5,0);        
hold on
plot(gcHistBinsPerMov{1},combGaussHist)
xlim([-2e-6 2e-6])
xlabel('Gaussian Curvature, 1/nm^2')
ylabel('Probability')
title(['Gaussian Curvature Distribution Averaged Over All Movies, n =' num2str(nMovies)])
saveThatShit(['Combined Gaussian Curvature Distribution'],combOutDir);

figure

combDiffHist = mean(meanDiffHist,1);
combDiffHistSTD = std(meanDiffHist,[],1);

plotTransparent(dpcHistBinsPerMov{1},combDiffHist,combDiffHistSTD,[0 0 1],.5,0);        
hold on
plot(dpcHistBinsPerMov{1},combDiffHist)
xlim([1e-4 1.9e-3])
xlabel('Difference of Principle Curvatures, 1/nm')
ylabel('Probability')
title(['Difference of Principle Curvatures Distribution Averaged Over All Movies, n =' num2str(nMovies)])
saveThatShit(['Combined Difference of Principle Curvatures Distribution'],combOutDir);

figure

combMAPCHist = mean(meanMAPCHist,1);
combMAPCHistSTD = std(meanMAPCHist,[],1);

plotTransparent(mapcHistBinsPerMov{1},combMAPCHist,combMAPCHistSTD,[0 0 1],.5,0);        
hold on
plot(mapcHistBinsPerMov{1},combMAPCHist)
xlim([0 5e-3])
xlabel('Maximum Absolute Principle Curvature, 1/nm')
ylabel('Probability')
title(['Maximum Absolute Principle Curvature Distribution Averaged Over All Movies, n =' num2str(nMovies)])
saveThatShit(['Combined Max Absolute Principle Curvature Distribution'],combOutDir);

%--------- All Frames All Movies Simple Stat Averages -------- %

allStat = horzcat(mapcMedPerMov{:});%Do it this way so cut-paste is easier. Man I am a shitty programmer.
currStat = 'Median of Maximum Absolute Curvature';
meanAll = mean(allStat);
stdAll = std(allStat);
figure
hist(allStat)
xlabel(currStat)
ylabel('Count, All Frames, All Movies')
title({[currStat ', All Frames, All Movies'],...
        ['Mean: ' num2str(meanAll) ' +/- ' num2str(stdAll) ' STD']});
saveThatShit([currStat ' all frames all movies'],combOutDir);

allStat = horzcat(mapcHistModePerMov{:});%Do it this way so cut-paste is easier. Man I am a shitty programmer.
currStat = 'Mode of Maximum Absolute Curvature';
meanAll = mean(allStat);
stdAll = std(allStat);
figure
hist(allStat)
xlabel(currStat)
ylabel('Count, All Frames, All Movies')
title({[currStat ', All Frames, All Movies'],...
        ['Mean: ' num2str(meanAll) ' +/- ' num2str(stdAll) ' STD']});
saveThatShit([currStat ' all frames all movies'],combOutDir);

allStat = horzcat(mapcMeanPerMov{:});%Do it this way so cut-paste is easier. Man I am a shitty programmer.
currStat = 'Mean of Maximum Absolute Curvature';
meanAll = mean(allStat);
stdAll = std(allStat);
figure
hist(allStat)
xlabel(currStat)
ylabel('Count, All Frames, All Movies')
title({[currStat ', All Frames, All Movies'],...
        ['Mean: ' num2str(meanAll) ' +/- ' num2str(stdAll) ' STD']});
saveThatShit([currStat ' all frames all movies'],combOutDir);

allStat = horzcat(mcMeanPerMov{:});%Do it this way so cut-paste is easier. Man I am a shitty programmer.
currStat = 'Mean of Mean Curvature';
meanAll = mean(allStat);
stdAll = std(allStat);
figure
hist(allStat)
xlabel(currStat)
ylabel('Count, All Frames, All Movies')
title({[currStat ', All Frames, All Movies'],...
        ['Mean: ' num2str(meanAll) ' +/- ' num2str(stdAll) ' STD']});
saveThatShit([currStat ' all frames all movies'],combOutDir);

allStat = horzcat(mcMedPerMov{:});%Do it this way so cut-paste is easier. Man I am a shitty programmer.
currStat = 'Median of Mean Curvature';
meanAll = mean(allStat);
stdAll = std(allStat);
figure
hist(allStat)
xlabel(currStat)
ylabel('Count, All Frames, All Movies')
title({[currStat ', All Frames, All Movies'],...
        ['Mean: ' num2str(meanAll) ' +/- ' num2str(stdAll) ' STD']});
saveThatShit([currStat ' all frames all movies'],combOutDir);


save([combOutDir filesep 'combined stats and hists all movies.mat'],'combMeanHist','combMeanHistSTD','histBinsPerMov',...
    'combGaussHist','combGaussHistSTD','gcHistBinsPerMov',...
    'combDiffHist','combDiffHistSTD','dpcHistBinsPerMov',...
    'combMAPCHist','combMAPCHistSTD','mapcHistBinsPerMov',...
    'mapcMeanPerMov','mapcMedPerMov','mapcVarPerMov',...
    'mcMeanPerMov','mcMedPerMov','mcVarPerMov',...
    'nCurvPtsPerMov')



