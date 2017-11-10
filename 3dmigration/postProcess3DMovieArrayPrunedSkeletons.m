function postProcess3DMovieArrayPrunedSkeletons(MA,varargin)
%POSTPROCESS3DMOVIEARRAYPRUNEDSKELETONS calculates various statistics regarding the pruned skeletons of the input movie array 
% 
% 
% 
% 
% 
% 
%TEMP - ADD SOME FUCKING DOCUMENTATION!!! 
% 
%TEMP - separate single-cell analysis to it's own function? Or at least
%save to each movies directory so other functions can use it!!!
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

perMovDirName = 'pruned skeleton post processing';%Directory for saving individual move results in movie output directory

outFileName = 'pruned skeleton post processing.mat';%For saving results that are calculated here
outVars = {};

%Print options for saving figures to disk
pOpt = {'-r300',...% dpi = 300
        '-depsc2'};% use eps format
pOptTIFF = {'-r100','-dtiff'};%150 dpi for TIF format since this is usally just used for emailing to Bob!
   

iProcChan = 1;

nBins2D = 50; %Number of bins for 2D histogram/pdfs
pct2D = 67.5; %Percentile to outline in 2D hist plots. This is ~ 1 STD   


%% ------------- Input ------------ %%


%Parse all the inputs
ip = inputParser;
ip.FunctionName = mfilename;
ip.addRequired('MA',@(x)(isa(x,'MovieData3D')));
ip.addParamValue('ChannelIndex',1,@(x)(numel(x) == 1 && isposint(x)));
ip.addParamValue('TipRadiusRange',[500 Inf],@(x)(numel(x) == 2 && diff(x) > 0));%Tip radius range to include in the simplifying "branch count" figure. This weeds out filopoida and also false-positives which are almost universally small.
ip.addParamValue('BatchMode',false,(@(x)(numel(x)==1)));
ip.addParamValue('OutputDirectory',[],@ischar);%Directory will be created if it doesn't exist
ip.parse(MA,varargin{:});

p = ip.Results;

if isempty(p.OutputDirectory)
    p.OutputDirectory = uigetdir(pwd,'Select a directory to store the results:');
    if p.OutputDirectory == 0
        error('You must select an output directory!')
    end
end
%Just increment if folder exists
mkClrDir(p.OutputDirectory);

%% ------------ Init --------------- %%


nMovies = numel(MA);

nFramesPerMov = nan(nMovies,1);

if ~exist(p.OutputDirectory,'dir')
    mkdir(p.OutputDirectory)
end

if ~p.BatchMode
    wtBar = waitbar(0,'Please wait, performing post-processing on each movie...');
end        


histBins = [1 3 4 5 6 7];%Skip two - think about it.
nBins = numel(histBins);
vdHists = cell(nMovies,1);%Vertex degree histograms for each movie

nAngBins = 10;%Number of bins for angle histograms

meanHist = nan(nMovies,nBins);
branchDir = cell(nMovies,1);
branchAng = cell(nMovies,1);

tipPathLen = cell(nMovies,1);
tipPathN = cell(nMovies,1);


phiBins = linspace(-pi/2,pi/2,nAngBins);
thetaBins = linspace(-pi,pi,nAngBins);


allTipPathLen = cell(nMovies,1);
allTipPathN = cell(nMovies,1);

%Fix these so we have common value across movies/cond
tipPathBinSz = 5;%in microns
tipRadBinSz = 250;%in nanometers

hasMG = false(nMovies,1);

branchRad = cell(nMovies,1);
branchDists = cell(nMovies,1);
hasDist = cell(nMovies,1);

%% ---------------- Per-Movie Processing ---------------- %%


for iMov = 1:nMovies
    
    iSkProc = MA(iMov).getProcessIndex('SkeletonPruningProcess',1,~p.BatchMode);
    iMgProc = MA(iMov).getProcessIndex('MaskGeometry3DProcess',1,~p.BatchMode);
    if ~isempty(iMgProc) && MA(iMov).processes_{iMgProc}.checkChannelOutput(iProcChan)
            hasMG(iMov) = true;
    end
    
    if ~isempty(iSkProc) && MA(iMov).processes_{iSkProc}.checkChannelOutput(p.ChannelIndex);
    
        if isempty(MA(iMov).eventTimes_)
            nFramesPerMov(iMov) = MA(iMov).nFrames_;
        else
            %TEMP - this doesn't currently support a start frame > 1!!!
            nFramesPerMov(iMov) = MA(iMov).eventTimes_(2);
        end
        
        vdHists{iMov} = nan(nFramesPerMov(iMov),nBins);
        
        currOutDir = [MA(iMov).outputDirectory_ filesep perMovDirName];
        mkClrDir(currOutDir);
        
        branchDir{iMov} = cell(nFramesPerMov(iMov),1);
        branchAng{iMov} = cell(nFramesPerMov(iMov),1);
        tipPathLen{iMov} = cell(nFramesPerMov(iMov),1);
        tipPathN{iMov} = cell(nFramesPerMov(iMov),1);
        branchRad{iMov} = cell(nFramesPerMov(iMov),1);
        branchDists{iMov} = cell(nFramesPerMov(iMov),1);
        hasDist{iMov} = cell(nFramesPerMov(iMov),1);
        
        for iFrame = 1:nFramesPerMov(iMov);
        
            currSkel = MA(iMov).processes_{iSkProc}.loadChannelOutput(p.ChannelIndex,iFrame);
            [iTipVert{iMov}{iFrame},iTipsEdge{iMov}{iFrame}] = findTips(currSkel.edges,size(currSkel.vertices,1));
            if hasMG(iMov)
                currMaskProp = MA(iMov).processes_{iMgProc}.loadChannelOutput(iProcChan,iFrame);
                branchRad{iMov}{iFrame} = branchRadii(currSkel,currMaskProp);
                branchRad{iMov}{iFrame} = cellfun(@(x)(x .* MA(iMov).pixelSize_),branchRad{iMov}{iFrame},'Unif',0);%Convert to nm
                
                meanBranchRad{iMov}{iFrame} = cellfun(@(x)(mean(x)),branchRad{iMov}{iFrame});
                tipRad{iMov}{iFrame} = branchRad{iMov}{iFrame}(iTipsEdge{iMov}{iFrame});
                meanTipRad{iMov}{iFrame} = meanBranchRad{iMov}{iFrame}(iTipsEdge{iMov}{iFrame});
            end
            
            % ---- Simplified "branch count" Analysis ---- %
            
            %Get the simplified thresholded branch count which excludes
            %small structures which are presumably filopodia or false
            %positives
            nTipsPerFrameThresh{iMov}(iFrame) = numel(meanTipRad{iMov}{iFrame} >= p.TipRadiusRange(1) & meanTipRad{iMov}{iFrame} <= p.TipRadiusRange(2));            
            nTipsPerFrame{iMov}(iFrame) = numel(iTipVert{iMov}{iFrame});
            
            % ----- Vertex Degree Analysis ----- %
            
            %skelStats = skelGraphStats(currSkel.vertices,currSkel.edges,currSkel.edgePaths);
            skelStats.vertexDegree = graphVertDegree(currSkel.edges,size(currSkel.vertices,1));                         
            vdHists{iMov}(iFrame,:) = histc(skelStats.vertexDegree(skelStats.vertexDegree > 0),histBins);
            vertDegreeMean{iMov}(iFrame) = mean(skelStats.vertexDegree);
            vertDegreeMedian{iMov}(iFrame) = median(skelStats.vertexDegree);
            vertBranchPointMean{iMov}(iFrame) = mean(skelStats.vertexDegree(skelStats.vertexDegree>2));
            vertBranchPointTotal{iMov}(iFrame) = sum(skelStats.vertexDegree(skelStats.vertexDegree>2));
            
            % ----- Branch Direction Analysis ----- %
            
            [branchDir{iMov}{iFrame},branchAng{iMov}{iFrame}] = calcBranchTipDirections(currSkel.edges,currSkel.edgePaths,size(currSkel.vertices,1));                                                                                                
            %And get the mean direction in each frame
            meanBranchDir{iMov}(iFrame,:) = nanmean(branchDir{iMov}{iFrame},1);
            [meanBranchAng{iMov}(iFrame,1) meanBranchAng{iMov}(iFrame,2) ~] = cart2sph(meanBranchDir{iMov}(iFrame,1),meanBranchDir{iMov}(iFrame,2),meanBranchDir{iMov}(iFrame,3));
            
           
            % ----- Tip Path Length Analysis ---- %
                        
            [tmp1,tmp2] = analyzeSkeletonTipPaths(currSkel.vertices,currSkel.edges,currSkel.edgePaths,currSkel.edgeLabels);
            %Convert to microns, workaround for occasional errors in skeleton
            %graphs
            tipPathLen{iMov}{iFrame} = tmp1(~isnan(tmp1) & ~isinf(tmp1)) .* MA(iMov).pixelSize_ / 1e3;
            tipPathN{iMov}{iFrame} = cellfun(@numel,tmp2(~isnan(tmp1) & ~isinf(tmp1)));                        
            
            %----- Distance From Centermost Point Analysis ---- %
            branchDists{iMov}{iFrame} = analyzeSkeletonDistanceFromPoint(currSkel.vertices,...
                                                                         currSkel.edges,...
                                                                         currSkel.edgePaths,...
                                                                         mean(currMaskProp.CenterMost(:,[2 1 3]),1));%Centermost is in XYZ, and may not be unique so we take mean                                                                                 
            %Longest edges in loops currently have undefined distance..
            %keep track of these
            hasDist{iMov}{iFrame} = ~cellfun(@isempty,branchDists{iMov}{iFrame});
            
            %Convert to nm
            branchDists{iMov}{iFrame}(hasDist{iMov}{iFrame}) = cellfun(@(x)(x .* MA(iMov).pixelSize_),branchDists{iMov}{iFrame}(hasDist{iMov}{iFrame}),'Unif',0);
        end        
        
        tData{iMov} = 0:MA(iMov).timeInterval_:(MA(iMov).timeInterval_*(nFramesPerMov(iMov)-1));
        
        % ----- All-frame mean calcs per movie --------- %
        
        %Branch-point degree
        meanHist(iMov,:) = mean(vdHists{iMov},1);
        
        %Branch direction
        
        %Get average directions. We can't average the angles themselves
        %because they're circular (0==2pi), so we average the vectors and
        %then get the angle
        
        %NOTE: This is not what we want - if the x/y are random and average
        %out, then any small z component bias gives an exagerrated
        %angle!!!!!!!! Need better way to average, or get mask body main axis
        %and compare to that!!
        
        allBranchDir{iMov} = vertcat(branchDir{iMov}{:});
        avgDir(iMov,:) = nanmean(allBranchDir{iMov},1);
        [avgAng(iMov,1),avgAng(iMov,2)] = cart2sph(avgDir(iMov,1),avgDir(iMov,2),avgDir(iMov,3));
        allBranchAng{iMov} = vertcat(branchAng{iMov}{:});
        branchAngHist{iMov}(:,2) = histc(vertcat(allBranchAng{iMov}(:,2)),phiBins);
        branchAngHist{iMov}(:,1) = histc(vertcat(allBranchAng{iMov}(:,1)),thetaBins);
        
        %Branch-tip paths
        allTipPathLen{iMov} = vertcat(tipPathLen{iMov}{:});
        allTipPathN{iMov} = vertcat(tipPathN{iMov}{:});
        allTipLenBins{iMov} = 0:tipPathBinSz:max(allTipPathLen{iMov});
        tipLenHist{iMov} = histc(allTipPathLen{iMov},allTipLenBins{iMov});
        allTipNBins{iMov} = 0:max(allTipPathN{iMov});
        tipNHist{iMov} = histc(allTipPathN{iMov},allTipNBins{iMov});        
        
        %Branch & Tip Radii
        allMeanTipRad{iMov} = vertcat(meanTipRad{iMov}{:});
        tipRadBins{iMov} = 0:tipRadBinSz:max(allMeanTipRad{iMov});
        
        %Thresholded tip count
        nTipsMeanPerMovThresh(iMov) = mean(nTipsPerFrameThresh{iMov});
        nTipsSTDPerMovThresh(iMov) = std(nTipsPerFrameThresh{iMov});        
                
        % Per-Movie Figures
        
        % ---- Vertex Degree Figure ----- %
        
        %TEMP - the x limits on this and the rest of the figures need to be
        %fixed! The auto-selection fucks it up most of the time...
        
        if p.BatchMode
            vdFigPM(iMov) = figure('Visible','off');
        else
            vdFigPM(iMov) = figure;
        end
        subplot(1,2,1);
        bar(meanHist(iMov,:))
        set(gca,'XTickLabel',histBins)
        xlabel('Branch Point Degree')
        ylabel('Average Count Per Frame')
        title('Branch point degree distribution, average per-frame')
        subplot(1,2,2);
        set(vdFigPM(iMov),'DefaultAxesColorOrder',jet(nFramesPerMov(iMov)));
        plot(vdHists{iMov}')
        set(gca,'XTickLabel',histBins)
        xlabel('Branch Point Degree')
        ylabel('Count in Each Frame')                
        figName = [currOutDir filesep 'Branch Degree Distribution'];
        print(pOpt{:},[figName '.eps']);
        print(pOptTIFF{:},[figName '.tif']);
        hgsave([figName '.fig'])
        
        if p.BatchMode
            vdotFigPM(iMov) = figure('Visible','off');
        else
            vdotFigPM(iMov) = figure;
        end
        
        plot(tData{iMov},vertDegreeMean{iMov})
        xlabel('Time, Seconds')
        ylabel('Mean Vertex Degree')
        ylim(ylim .* [.95 1.05])
        saveThatShit('mean vertex degree over time',currOutDir)
        
        if p.BatchMode
            vdotFigPM(iMov) = figure('Visible','off');
        else
            vdotFigPM(iMov) = figure;
        end
        
        plot(tData{iMov},vertBranchPointMean{iMov})
        xlabel('Time, Seconds')
        ylabel('Mean Branch Point Degree')
        ylim(ylim .* [.95 1.05])
        saveThatShit('mean branch point degree over time',currOutDir)
        
        if p.BatchMode
            vdotFigPM(iMov) = figure('Visible','off');
        else
            vdotFigPM(iMov) = figure;
        end
        
        plot(tData{iMov},vertBranchPointTotal{iMov})
        xlabel('Time, Seconds')
        ylabel('Total Branch Point Degree')
        ylim(ylim .* [.95 1.05])
        saveThatShit('total branch point degree over time',currOutDir)
                                            
        
        % -------- Branch Direction Figure ------ %
        
        if p.BatchMode
            bdFigPM(iMov) = figure('Visible','off');
        else
            bdFigPM(iMov) = figure;
        end
                
        frameCols = jet(nFramesPerMov(iMov));               
        for j = 1:nFramesPerMov(iMov)      
                       
            thetaCts = histc(branchAng{iMov}{j}(:,1),thetaBins);
            phiCts = histc(branchAng{iMov}{j}(:,2),phiBins);
  
            
            subplot(4,1,1)
            if j == 1
                hold on
                xlabel('Phi, Radians')
                xlim([-pi/2 pi/2])
                ylabel('Count')
                title('Per-Frame Distribution, red=last frame, blue = first')
            end            
            plot(phiBins,phiCts,'color',frameCols(j,:));
            subplot(4,1,2)
            if j == 1
                hold on
                xlabel('Theta, Radians')
                ylabel('Count')
                xlim([-pi pi])
            end            
            plot(thetaBins,thetaCts,'color',frameCols(j,:));                                        
        end        
        
        subplot(4,1,3)
        hold on
        title('Whole-Movie Distributions')
        xlabel('Phi, Radians')
        ylabel('Count')
        xlim([-pi/2 pi/2])        
        bar(phiBins,branchAngHist{iMov}(:,2))        
        plot(avgAng(iMov,2) * ones(1,2),ylim,'--r')
        
        subplot(4,1,4)
        hold on        
        xlabel('Theta, Radians')
        ylabel('Count')       
        xlim([-pi pi])        
        bar(thetaBins,branchAngHist{iMov}(:,1))        
        plot(avgAng(iMov) * ones(1,2),ylim,'--r')                        
        
        figName = [currOutDir filesep 'Branch Angle Distribution'];
        print(pOpt{:},[figName '.eps']);
        print(pOptTIFF{:},[figName '.tif']);
        hgsave([figName '.fig'])
        
        if p.BatchMode
            baFig(iMov) = figure('Visible','off');
        else
            baFig(iMov) = figure;
        end
        h = rose(allBranchAng{iMov}(:,1),thetaBins);
        title('Absolute branch angle rose histogram over all frames')
        figName = [currOutDir filesep 'Absolute branch angle roseplot'];
        mfFigureExport(baFig(iMov),figName)
        
        if p.BatchMode
            batFig(iMov) = figure('Visible','off');
        else
            batFig(iMov) = figure;
        end
        hold on
        frameCols = jet(nFramesPerMov(iMov));
        for j = 1:nFramesPerMov(iMov)
            nBrCurr = size(branchDir{iMov}{j},1);
            quiver3(zeros(nBrCurr,1),zeros(nBrCurr,1),zeros(nBrCurr,1),branchDir{iMov}{j}(:,1),branchDir{iMov}{j}(:,2),branchDir{iMov}{j}(:,3),0,'color',frameCols(j,:));
        end
        title({'All branch directions over time',...
            'Color indicates frame time in seconds'})
        if nFramesPerMov(iMov) > 1
            colorbar
            caxis([min(tData{iMov}) max(tData{iMov})])
        end
        view(3)
        figName = [currOutDir filesep 'All branch directions over time 3D plot'];
        mfFigureExport(batFig(iMov),figName) 
        
        if p.BatchMode
            abdtFig(iMov) = figure('Visible','off');
        else
            abdtFig(iMov) = figure;
        end
        hold on
        for j = 1:nFramesPerMov(iMov)
        
            quiver3(0,0,0,meanBranchDir{iMov}(j,1),meanBranchDir{iMov}(j,2),meanBranchDir{iMov}(j,3),0,'color',frameCols(j,:))
            
        end
        view(3)
        if nFramesPerMov(iMov) > 1
            colorbar
            caxis([min(tData{iMov}) max(tData{iMov})])
        end
        title({'Mean branch direction over time',...
               'color indicates time, seconds'})
        figName = [currOutDir filesep 'Mean branch direction over time 3D plot'];
        mfFigureExport(abdtFig(iMov),figName)     
        
        if p.BatchMode
            abatFig(iMov) = figure('Visible','off');
        else
            abatFig(iMov) = figure;
        end
        plot(tData{iMov},unwrap(meanBranchAng{iMov}(:,1)))        
        xlabel('Time, Seconds')
        ylabel('Unwrapped mean branch theta angle, radians')
        title('Unwrapped mean branch theta over time')
        figName = [currOutDir filesep 'Mean branch angle unwrapped over time 3D plot'];
        mfFigureExport(abatFig(iMov),figName)     
        
        if p.BatchMode
            abatFig(iMov) = figure('Visible','off');
        else
            abatFig(iMov) = figure;
        end
        plot(tData{iMov},meanBranchAng{iMov}(:,1))        
        xlabel('Time, Seconds')
        ylabel('Mean branch theta angle, radians')
        title('Mean branch theta over time')
        hold on
        plot(xlim,[pi pi],'--k')
        plot(xlim,[-pi -pi],'--k')
        figName = [currOutDir filesep 'Mean branch angle over time 3D plot'];
        mfFigureExport(abatFig(iMov),figName)     
        
        
        
        % ----- Tip Path Figure ---- %
        if p.BatchMode
            tpFigPM(iMov) = figure('Visible','off');
        else
            tpFigPM(iMov) = figure;
        end
        subplot(2,1,1)
        bar(allTipLenBins{iMov},tipLenHist{iMov})
        hold on
        xlabel('Path Length, Tip To cell Body (microns)')
        ylabel('Count')
        title('Tip-To-Cell-Body Path Length Distribution, All frames')
        
        subplot(2,1,2)
        bar(allTipNBins{iMov},tipNHist{iMov})
        hold on
        xlabel('# Of vertices on path')
        ylabel('Count')
        title('Tip-To-Cell-Body Path Complexity Distribution, All frames')
        
        figName = [currOutDir filesep 'Tip Path Distributions'];
        print(pOpt{:},[figName '.eps']);
        print(pOptTIFF{:},[figName '.tif']);
        hgsave([figName '.fig'])
        
        % ----- Radius vs. Distance from Centermost Figure ---- %
        if p.BatchMode
            drFigPM(iMov) = figure('Visible','off');
        else
            drFigPM(iMov) = figure;                       
        end
        subplot(2,1,1)
        hold on
        currAllDists = vertcat(branchDists{iMov}{:});
        currHasDists = vertcat(hasDist{iMov}{:});
        currAllDists = vertcat(currAllDists{:});        
        currAllRad = vertcat(branchRad{iMov}{:});
        currAllRad = vertcat(currAllRad{currHasDists});
        [N,C] = hist3([currAllDists currAllRad],[20 20]);
        imagesc(C{1},C{2},log(N')),axis xy        
        xlim([min(C{1}) max(C{1})])
        ylim([min(C{2}) max(C{2})])
        title({'Distance along skeleton from centermost point',...
                'vs. radius, Log10 of histogram',...
                'all frames'});
        xlabel('Distance, nm')
        ylabel('Radius, nm')        
        
        subplot(2,1,2)
        hold on
        for j = 1:nFramesPerMov(iMov)                        
            plot(vertcat(branchDists{iMov}{j}{hasDist{iMov}{j}}),vertcat(branchRad{iMov}{j}{hasDist{iMov}{j}}),'.','color',frameCols(j,:))
        end
        title({'Distance along skeleton from centermost point',...
                'vs. radius, per frame',...
                'Blue=1st frame, Red=Last'});
        xlabel('Distance, nm')
        ylabel('Radius, nm')
        
        
        figName = [currOutDir filesep 'radius vs distance from centermost point'];
        print(pOpt{:},[figName '.eps']);
        print(pOptTIFF{:},[figName '.tif']);
        hgsave([figName '.fig'])
        
        
        % ----- Branch tip radius histogram ---- %
        if p.BatchMode
            trFigPM(iMov) = figure('Visible','off');
        else
            trFigPM(iMov) = figure;                       
        end
        
        tipRadHistMeanPerFramePerMov{iMov} = histc(allMeanTipRad{iMov},tipRadBins{iMov});
        tipRadHistMeanPerFramePerMov{iMov} = tipRadHistMeanPerFramePerMov{iMov} ./ nFramesPerMov(iMov);
        bar(tipRadBins{iMov},tipRadHistMeanPerFramePerMov{iMov})
        xlabel('Tip Radius, nm')
        ylabel('Mean Count, Per Frame')
        title('Mean Tip Radius Histogram, Per Frame')
        
        figName = [currOutDir filesep 'tip radius histogram'];
        print(pOpt{:},[figName '.eps']);
        print(pOptTIFF{:},[figName '.tif']);
        hgsave([figName '.fig'])
        
        % ----- Thresholded tip count figure ----- %
        
        if nFramesPerMov(iMov) > 1 
            if p.BatchMode
                ttcFigPM(iMov) = figure('Visible','off');
            else
                ttcFigPM(iMov) = figure;                       
            end
            hold on
            title({'Thresholded Tip Count vs. Time',['Tip Radii ' num2str(p.TipRadiusRange(1)) ' < x < ' num2str(p.TipRadiusRange(2)) ' nm'],...
                ['Mean, all frames: ' num2str(nTipsMeanPerMovThresh(iMov)) ', STD = ' num2str(nTipsSTDPerMovThresh(iMov))]})
            plot(0:MA(iMov).timeInterval_:(MA(iMov).timeInterval_*(nFramesPerMov(iMov)-1)),nTipsPerFrameThresh{iMov},'.-')
            xlabel('Time, Seconds')
            ylabel('Thresholded Tip Count')

            figName = [currOutDir filesep 'tip count thresholded vs time'];
            print(pOpt{:},[figName '.eps']);
            print(pOptTIFF{:},[figName '.tif']);
            hgsave([figName '.fig'])
        end
        %TEMP - Save results to individual movie output folders also!!!!
        
    
    else
        warning('MIGRATION3D:ppSkelGraph:noPruning',...
            ['Movie ' num2str(iMov) ' does not have valid skeleton pruning - not analyzing!']);
    end
    
    if ~p.BatchMode    
        if ishandle(wtBar)
            waitbar(iMov/nMovies,wtBar)
        end
        close all
    else
        close all
    end
        
end


outVars = [outVars 'nTipsPerFrame'];

if nnz(hasMG) ~=nMovies
    warning('Some movies did not have valid mask geometry analysis!!!');
end
if ishandle(wtBar)
    close(wtBar);
end

% ------ Combined (All-Movie) Branch-Point-Degree Data ------- %

combMean = nanmean(meanHist,1);
if nMovies > 1
    combSTD = nanstd(meanHist,1);
else
    combSTD = zeros(1,numel(combMean));
end

if p.BatchMode
    cmFig(iMov) = figure('Visible','off');
else
    cmFig(iMov) = figure;
end
subplot(2,1,1)
bar(combMean)
hold on
errorbar(combMean,combSTD * 1.96 ./ sqrt(nMovies),'.r')
set(gca,'XTickLabel',histBins)
xlabel('Branch Point Degree')
ylabel('Average Count Per Frame')
title('Branch point degree distribution, average per-frame and ~95% C.I.')
subplot(2,1,2)
combProbDist = nanmean(meanHist ./ repmat(sum(meanHist,2),[1 numel(histBins)]),1);
combProbSTD = nanstd(meanHist ./ repmat(sum(meanHist,2),[1 numel(histBins)]),[],1);
bar(combProbDist)
hold on
errorbar(combProbDist,combProbSTD * 1.96 ./ sqrt(nMovies),'.r')
set(gca,'XTickLabel',histBins)
xlabel('Branch Point Degree')
ylabel('Average Probability')
title({'Branch point degree probability distribution, all movies & ~95% C.I.',['n=' num2str(nMovies)]})


figName = [p.OutputDirectory filesep 'Combined Branch Degree Distribution'];
print(pOpt{:},[figName '.eps']);
print(pOptTIFF{:},[figName '.tif']);
hgsave([figName '.fig'])


% ------------- Combined Branch-Direction Data ------------ %


combAngData = vertcat(allBranchAng{:});
combPhiHist = histc(combAngData(:,2),phiBins);
meanPhi = nanmean(combAngData(:,2));
stdPhi = nanstd(combAngData(:,2));
figure
subplot(2,1,1)
bar(phiBins,combPhiHist)
xlim([-pi/2 pi/2])
xlabel('Phi, Radians')
ylabel('Count')
title({'Unaligned, combined Phi Histogram, all Movies',['Combined Phi STD: ' num2str(stdPhi)],['n=' num2str(nMovies)]})

combThetaHist = histc(combAngData(:,1),thetaBins);
meanTheta = nanmean(combAngData(:,1));
stdTheta = nanstd(combAngData(:,1));
subplot(2,1,2);
bar(thetaBins,combThetaHist)
xlim([-pi pi])
xlabel('Theta, Radians')
ylabel('Count')
title({'Unaligned, combined Theta Histogram, all Movies',['Combined Theta STD: ' num2str(stdTheta)],['n=' num2str(nMovies)]})

figName = [p.OutputDirectory filesep 'Combined Branch Angle Distributions'];
print(pOpt{:},[figName '.eps']);
print(pOptTIFF{:},[figName '.tif']);
hgsave([figName '.fig'])




%% --------------------- Combined Tip-Path Data ------------------- %


% --- PDFs ---- %

combPathLenBins = 0:tipPathBinSz:max(cellfun(@max,allTipLenBins));
combPathPDF = zeros(nMovies,numel(combPathLenBins));

combPathNBins = 0:max(cellfun(@max,allTipNBins));
combPathNPDF = zeros(nMovies,numel(combPathNBins));


for iMov = 1:nMovies
    
    nCur = numel(tipLenHist{iMov});    
    combPathPDF(iMov,1:nCur) = combPathPDF(iMov,1:nCur) + (tipLenHist{iMov} ./ sum(tipLenHist{iMov}))';

    nCur=  numel(tipNHist{iMov});
    combPathNPDF(iMov,1:nCur) = combPathNPDF(iMov,1:nCur) + (tipNHist{iMov} ./ sum(tipNHist{iMov}))';
    
end

figure
subplot(2,1,1)
bar(combPathLenBins,mean(combPathPDF,1))
hold on
errorbar(combPathLenBins,mean(combPathPDF,1),std(combPathPDF,[],1)*1.96 / sqrt(nMovies),'.r')
maxBin = max(combPathLenBins);
xlim([-tipPathBinSz/2 maxBin + tipPathBinSz/2])
yMax = max(ylim);
ylim([0 yMax])
title({'Combined Tip-To-Body Path-Length Probability Distribution, All Movies',['n=' num2str(nMovies)]})
xlabel('Path Length, microns')
ylabel('Average Probability')

subplot(2,1,2)
bar(combPathNBins,mean(combPathNPDF,1))
hold on
errorbar(combPathNBins,mean(combPathNPDF,1),std(combPathNPDF,[],1) * 1.96 / sqrt(nMovies),'.r')
maxBin = max(combPathNBins);
xlim([1.5 maxBin + 1/2])
yMax = max(ylim);
ylim([0 yMax])
title('Combined Tip-To-Body Path Complexity Probability Distribution, All Movies')
xlabel('Path Complexity, # of Vertices')
ylabel('Average Probability')

figName = [p.OutputDirectory filesep 'Combined Tip Path Distributions'];
print(pOpt{:},[figName '.eps']);
print(pOptTIFF{:},[figName '.tif']);
hgsave([figName '.fig'])


outVars = [outVars 'combPathLenBins','combPathPDF','combPathNBins','combPathNPDF','combPathHist','combPathNHist'];

%% --- Avg Per Mov Path Length Histograms --- %%

combPathLenBins = 0:tipPathBinSz:max(cellfun(@max,allTipLenBins));
combPathHist = zeros(nMovies,numel(combPathLenBins));

combPathNBins = 0:max(cellfun(@max,allTipNBins));
combPathNHist = zeros(nMovies,numel(combPathNBins));


for iMov = 1:nMovies
    
    nCur = numel(tipLenHist{iMov});    
    combPathHist(iMov,1:nCur) = combPathHist(iMov,1:nCur) + (tipLenHist{iMov} ./ nFramesPerMov(iMov))';

    nCur=  numel(tipNHist{iMov});
    combPathNHist(iMov,1:nCur) = combPathNHist(iMov,1:nCur) + (tipNHist{iMov} ./ nFramesPerMov(iMov))';
    
end

figure
subplot(2,1,1)
bar(combPathLenBins,mean(combPathHist,1))
hold on
errorbar(combPathLenBins,mean(combPathHist,1),std(combPathHist,[],1)*1.96 / sqrt(nMovies),'.r')
title({'Combined Tip-To-Body Path-Length Mean Per-Frame Histogram, All Movies',['n=' num2str(nMovies)]})
xlabel('Path Length, microns')
ylabel('Average Count Per-Frame')

subplot(2,1,2)
bar(combPathNBins,mean(combPathNHist,1))
hold on
errorbar(combPathNBins,mean(combPathNHist,1),std(combPathNHist,[],1) * 1.96 / sqrt(nMovies),'.r')
title('Combined Tip-To-Body Path Complexity Mean Per-Frame Histogram, All Movies')
xlabel('Path Complexity, # of Vertices')
ylabel('Average Count Per-Frame')

figName = [p.OutputDirectory filesep 'Combined Tip Path Histograms'];
print(pOpt{:},[figName '.eps']);
print(pOptTIFF{:},[figName '.tif']);
hgsave([figName '.fig'])


outVars = [outVars 'combPathLenBins','combPathHist','combPathNBins','combPathNHist','combPathHist','combPathNHist'];


% --------------------- Combined Tip Radius Data ------------------- %

combTipRadBins =  0:tipRadBinSz:max(cellfun(@max,allMeanTipRad));
combTipRadHist = zeros(nMovies,numel(combTipRadBins));


for iMov = 1:nMovies
    
    nCur = numel(tipRadHistMeanPerFramePerMov{iMov});    
    combTipRadHist(iMov,1:nCur) = tipRadHistMeanPerFramePerMov{iMov};
    
end

figure
bar(combTipRadBins,mean(combTipRadHist,1))
hold on
errorbar(combTipRadBins,mean(combTipRadHist,1),std(combTipRadHist,[],1)*1.96 / sqrt(nMovies),'.r')
title({'Combined Tip Radius Mean Per-Frame Histogram, All Movies',['n=' num2str(nMovies)]})
xlabel('Tip Radius, nm')
ylabel('Average Count Per-Frame')

figName = [p.OutputDirectory filesep 'Combined Tip Radius Histogram'];
print(pOpt{:},[figName '.eps']);
print(pOptTIFF{:},[figName '.tif']);
hgsave([figName '.fig'])


outVars = [outVars 'combTipRadBins','combTipRadHist'];



%% -------------- Combined Radius Vs. Distance Data --------------- %


% ------ Radius Vs. Dist Hist ----- %

figure
hold on
allDists = cell(nMovies,1);
allHasDists = cell(nMovies,1);
allRad = cell(nMovies,1);


for iMov = 1:nMovies
    allDists{iMov} = vertcat(branchDists{iMov}{:});
    allHasDists{iMov} = vertcat(hasDist{iMov}{:});
    allDists{iMov} = vertcat(allDists{iMov}{:});        
    allRad{iMov} = vertcat(branchRad{iMov}{:});
    allRad{iMov} = vertcat(allRad{iMov}{allHasDists{iMov}});
end

allDists = vertcat(allDists{:});
allRad = vertcat(allRad{:});

histEdgesDist = linspace(min(allDists(:)),max(allDists(:))+eps,nBins2D+1);
histEdgesRad = linspace(min(allRad(:)),max(allRad(:))+eps,nBins2D+1);

[distRad2DHist,distRad2DBins] = hist3([allDists allRad],'Edges',{histEdgesDist histEdgesRad});

meanAllRadPerDist = nan(nBins2D,1);
rangeAllRadPerDist = nan(nBins2D,2);
for j = 1:nBins2D
    
    meanAllRadPerDist(j) = nanmean(allRad(allDists >= histEdgesDist(j) & allDists < histEdgesDist(j+1)));
    rangeAllRadPerDist(j,:) = prctile(allRad(allDists >= histEdgesDist(j) & allDists < histEdgesDist(j+1)),[(100-pct2D)/2  100-(100-pct2D)/2]);
    
end
                
imagesc(distRad2DBins{1},distRad2DBins{2},log(distRad2DHist')),axis xy
xlim([min(distRad2DBins{1}) max(distRad2DBins{1})])
ylim([min(distRad2DBins{2}) max(distRad2DBins{2})])
title({'Distance along skeleton from centermost point',...
        'vs. radius, Log10 of histogram',...
        ['all cells, n=' num2str(nMovies)]});
xlabel('Distance, nm')
ylabel('Radius, nm')        
colormap gray
hold on
plot(distRad2DBins{1}(1:end-1),meanAllRadPerDist,'r');
plot(distRad2DBins{1}(1:end-1),rangeAllRadPerDist(:,1),'--r');
legend('Mean',['Center ' num2str(pct2D) '%']);
plot(distRad2DBins{1}(1:end-1),rangeAllRadPerDist(:,2),'--r');

figName = [p.OutputDirectory filesep 'Combined Distance vs Radius Log Histogram'];
print(pOpt{:},[figName '.eps']);
print(pOptTIFF{:},[figName '.tif']);
hgsave([figName '.fig'])

% ------ Radius Vs. Dist PDF and Isovals ----- %
pPDF.xy = vertcat(distRad2DBins{:});
distRad2Dpdf = gkde2([allDists allRad],pPDF);

figure
hold on
imagesc(distRad2DBins{1},distRad2DBins{2},distRad2Dpdf.pdf)
xlim([min(distRad2DBins{1}) max(distRad2DBins{1})])
ylim([min(distRad2DBins{2}) max(distRad2DBins{2})])
title({'Distance along skeleton from centermost point',...
        'vs. radius, 2D PDF Estimate',...
        ['all cells, n=' num2str(nMovies)]});
xlabel('Distance, nm')
ylabel('Radius, nm')        
colormap gray
hold on
colorbar
%Numerically find an isovalue which contains 95% of the distribution
isoVal = .90;
tryIsoVals = linspace(min(distRad2Dpdf.pdf(:)),max(distRad2Dpdf.pdf(:)),1e3);
cumProbPerVal = arrayfun(@(x)(sum(distRad2Dpdf.pdf(distRad2Dpdf.pdf(:)>x))),tryIsoVals);
[minErr,iBestIso] = min(abs(sum(distRad2Dpdf.pdf(:)) * isoVal -  cumProbPerVal));
distRad2Dpdf.isoValPct = isoVal;
distRad2Dpdf.isoValProb = tryIsoVals(iBestIso);
distRad2Dpdf.isoCont = separateContours(contourc(distRad2DBins{1}(1:end-1),distRad2DBins{2}(1:end-1),distRad2Dpdf.pdf,distRad2Dpdf.isoValProb*ones(1,2)));
cellfun(@(x)(plot(x(1,:),x(2,:),'r')),distRad2Dpdf.isoCont);
legend([num2str(isoVal*100) '% Isocontour'])
outVars = [outVars 'distRad2Dpdf'];

figName = [p.OutputDirectory filesep 'Combined Distance vs Radius PDF'];
print(pOpt{:},[figName '.eps']);
print(pOptTIFF{:},[figName '.tif']);
hgsave([figName '.fig'])


%% ------ The over-simplifying "number of branches per cell" figure ------ %%


figure
hold on
hist(nTipsMeanPerMovThresh)
xlabel('Thresholded Branch Tip Number')
ylabel('Count, Mean Per Frame all Movies')
title({['Thresholded Mean Per-Frame Branch Tip Count, All Movies n=' num2str(nMovies)],...
        ['Tip Radii ' num2str(p.TipRadiusRange(1)) ' < x < ' num2str(p.TipRadiusRange(2)) ' nm'],...
        ['Combined mean = ' num2str(mean(nTipsMeanPerMovThresh)) ', STD = ' num2str(std(nTipsMeanPerMovThresh))]});    

outVars = [outVars 'nTipsMeanPerMovThresh','nTipsSTDPerMovThresh'];
    
figName = [p.OutputDirectory filesep 'Combined Tip Count Thresholded'];
print(pOpt{:},[figName '.eps']);
print(pOptTIFF{:},[figName '.tif']);
hgsave([figName '.fig'])

%% -------- Branch # Variability Per Cell vs. Mean Variability Accross Cells Figure ------ %%

currFig = figure;

%Get per-movie % coefficient of variation
cvPctPerMovnTipsThresh = nTipsSTDPerMovThresh ./ nTipsMeanPerMovThresh * 100;
meannTipsThreshAll = mean(nTipsMeanPerMovThresh);
stdnTipsThreshAll = std(nTipsMeanPerMovThresh);
cvPctAllMovnTipsThresh = stdnTipsThreshAll / meannTipsThreshAll * 100;

hist(cvPctPerMovnTipsThresh)
hold on
xlabel('% Coefficient of Variation (STD/Mean)')
ylabel('# Cells')
title({'CV%, Thresholded Branch Tip #, per cell average/STD over all frames',...
    ['n = ' num2str(nMovies) ' cells'],...
    ['Mean = ' num2str(mean(cvPctPerMovnTipsThresh))],...
    ['Population CV% of per-cell means = ' num2str(cvPctAllMovnTipsThresh) '%']});


figName = [p.OutputDirectory filesep 'cv percent per cell and population thresholded tip count'];
mfFigureExport(currFig,figName);

outVars = [outVars 'cvPctPerMovnTipsThresh','cvPctAllMovnTipsThresh','meannTipsThreshAll','stdnTipsThreshAll'];

%% -------- Branch # Variability over Time vs. Across Cells Figure ------ %%

currFig = figure;
hold on

movColors = rand(nMovies,3);

allTData = arrayfun(@(iMov)(0:MA(iMov).timeInterval_:(MA(iMov).timeInterval_*(nFramesPerMov(iMov)-1))),1:nMovies,'Unif',false);
maxT = max(horzcat(allTData{:}));

for iMov = 1:nMovies
    %Do the mean/std bands per movie first so they're in the background
    mP = nTipsMeanPerMovThresh(iMov) + nTipsSTDPerMovThresh(iMov);
    mM = nTipsMeanPerMovThresh(iMov) - nTipsSTDPerMovThresh(iMov);
    fill([0 maxT maxT 0],[mP mP mM mM],movColors(iMov,:),'FaceAlpha',.2,'EdgeColor','none');
    
end


for iMov = 1:nMovies
     
     
    plot(allTData{iMov},nTipsPerFrameThresh{iMov},'.-','Color',movColors(iMov,:))
            
    %nTipsMeanPerMovThresh(iMov) = mean(nTipsPerFrameThresh{iMov});
    %nTipsSTDPerMovThresh(iMov) = std(nTipsPerFrameThresh{iMov});            
    
end




xlabel('Time, Seconds')
ylabel('Thresholded Tip Count')    

title({'Thresholded Tip Count vs. Time',['Tip Radii ' num2str(p.TipRadiusRange(1)) ' < x < ' num2str(p.TipRadiusRange(2)) ' nm'],...
        'Bands show mean + / - STD'});        
mfFigureExport(currFig,[p.OutputDirectory filesep 'thresholded branch count all cell overlay with mean and std']);

outVars = [outVars 'allTData' 'nTipsPerFrameThresh' 'nTipsSTDPerMovThresh','p'];

%% ---------- Save Output ------ %%

if p.BatchMode
    close all
end


outFile = [p.OutputDirectory filesep outFileName];

outVars = [outVars {'distRad2DHist','distRad2DBins','meanAllRadPerDist','rangeAllRadPerDist','pct2D','histEdgesDist','histEdgesRad'}];

save(outFile,outVars{:});




