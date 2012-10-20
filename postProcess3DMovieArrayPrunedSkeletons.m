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
pct2D = 67.5; %Percentile to outline in 2D hist plots. This is ~ STD
    
%% ------------- Input ------------ %%


%Parse all the inputs
ip = inputParser;
ip.FunctionName = mfilename;
ip.addRequired('MA',@(x)(isa(x,'MovieData3D')));
ip.addParamValue('ChannelIndex',1,@(x)(numel(x) == 1 && isposint(x)));
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
tipPathBinSz = 5;
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
            if hasMG(iMov)
                currMaskProp = MA(iMov).processes_{iMgProc}.loadChannelOutput(iProcChan,iFrame);
                branchRad{iMov}{iFrame} = branchRadii(currSkel,currMaskProp);
                branchRad{iMov}{iFrame} = cellfun(@(x)(x .* MA(iMov).pixelSize_),branchRad{iMov}{iFrame},'Unif',0);%Convert to nm
            end
            
            % ----- Vertex Degree Analysis ----- %
            
            %skelStats = skelGraphStats(currSkel.vertices,currSkel.edges,currSkel.edgePaths);
            skelStats.vertexDegree = graphVertDegree(currSkel.edges,size(currSkel.vertices,1));             
            vdHists{iMov}(iFrame,:) = histc(skelStats.vertexDegree(skelStats.vertexDegree > 0),histBins);
            
            % ----- Branch Direction Analysis ----- %
            
            [branchDir{iMov}{iFrame},branchAng{iMov}{iFrame}] = calcBranchTipDirections(currSkel.edges,currSkel.edgePaths,size(currSkel.vertices,1));                                                                                                
            
            
            % ----- Tip Path Length Analysis ---- %
                        
            [tmp1,tmp2] = analyzeSkeletonTipPaths(currSkel.vertices,currSkel.edges,currSkel.edgePaths,currSkel.edgeLabels);
            %Convert to nm, workaround for occasional errors in skeleton
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
        
        
        
        %TEMP - Save results to individual movie output folders also!!!!
        
    
    else
        warning('MIGRATION3D:ppSkelGraph:noPruning',...
            ['Movie ' num2str(iMov) ' does not have valid skeleton pruning - not analyzing!']);
    end
    
    if ~p.BatchMode        
        waitbar(iMov/nMovies,wtBar)
    else
        close all
    end
        
end

if nnz(hasMG) ~=nMovies
    warning('Some movies did not have valid mask geometry analysis!!!');
end

close(wtBar);

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




% --------------------- Combined Tip-Path Data ------------------- %

combPathLenBins = 0:tipPathBinSz:max(cellfun(@max,allTipLenBins));
combPathHist = zeros(nMovies,numel(combPathLenBins));

combPathNBins = 0:max(cellfun(@max,allTipNBins));
combPathNHist = zeros(nMovies,numel(combPathNBins));


for iMov = 1:nMovies
    
    nCur = numel(tipLenHist{iMov});    
    combPathHist(iMov,1:nCur) = combPathHist(iMov,1:nCur) + (tipLenHist{iMov} ./ sum(tipLenHist{iMov}))';

    nCur=  numel(tipNHist{iMov});
    combPathNHist(iMov,1:nCur) = combPathNHist(iMov,1:nCur) + (tipNHist{iMov} ./ sum(tipNHist{iMov}))';
    
end

figure
subplot(2,1,1)
bar(combPathLenBins,mean(combPathHist,1))
hold on
errorbar(combPathLenBins,mean(combPathHist,1),std(combPathHist,[],1)*1.96 / sqrt(nMovies),'.r')
title({'Combined Tip-To-Body Path-Length Probability Distribution, All Movies',['n=' num2str(nMovies)]})
xlabel('Path Length, microns')
ylabel('Average Probability')

subplot(2,1,2)
bar(combPathNBins,mean(combPathNHist,1))
hold on
errorbar(combPathNBins,mean(combPathNHist,1),std(combPathNHist,[],1) * 1.96 / sqrt(nMovies),'.r')
title('Combined Tip-To-Body Path Complexity Probability Distribution, All Movies')
xlabel('Path Complexity, # of Vertices')
ylabel('Average Probability')

figName = [p.OutputDirectory filesep 'Combined Tip Path Distributions'];
print(pOpt{:},[figName '.eps']);
print(pOptTIFF{:},[figName '.tif']);
hgsave([figName '.fig'])


outVars = [outVars{:} {'combPathLenBins','combPathHist','combPathNBins','combPathNHist','combPathHist','combPathNHist'}];

%% -------------- Combined Radius Vs. Distance Data --------------- %

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

figName = [p.OutputDirectory filesep 'Combined Distance vs Radius Distribution'];
print(pOpt{:},[figName '.eps']);
print(pOptTIFF{:},[figName '.tif']);
hgsave([figName '.fig'])


if p.BatchMode
    close all
end

%% ---------- Save Output ------ %%

outFile = [p.OutputDirectory filesep outFileName];

outVars = [outVars{:} {'distRad2DHist','distRad2DBins','meanAllRadPerDist','rangeAllRadPerDist','pct2D','histEdgesDist','histEdgesRad'}];

save(outFile,outVars{:});




