function postProcess3DMovieArrayPrunedSkeletons(MA,varargin)
%POSTPROCESS3DMOVIEARRAYPRUNEDSKELETONS calculates various statistics regarding the pruned skeletons of the input movie array 
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

perMovDirName = 'pruned skeleton post processing';%Directory for saving individual move results in movie output directory

%Print options for saving figures to disk
pOpt = {'-r300',...% dpi = 300
        '-depsc2'};% use eps format

%% ------------- Input ------------ %%


%Parse all the inputs
ip = inputParser;
ip.FunctionName = mfilename;
ip.addRequired('MA',@(x)(isa(x,'MovieData3D')));
ip.addParamValue('ChannelIndex',1,@(x)(numel(x) == 1 && isposint(x)));
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

nFramesPerMov = nan(nMovies,1);



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

%% ---------------- Per-Movie Processing ---------------- %%


for iMov = 1:nMovies
    
    iSkProc = MA(iMov).getProcessIndex('SkeletonPruningProcess',1,~p.BatchMode);
    
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
        
        for iFrame = 1:nFramesPerMov(iMov);
        
            currSkel = MA(iMov).processes_{iSkProc}.loadChannelOutput(p.ChannelIndex,iFrame);
            
            
            % ----- Vertex Degree Analysis ----- %
            
            %skelStats = skelGraphStats(currSkel.vertices,currSkel.edges,currSkel.edgePaths);
            skelStats.vertexDegree = graphVertDegree(currSkel.edges,size(currSkel.vertices,1));             
            vdHists{iMov}(iFrame,:) = histc(skelStats.vertexDegree(skelStats.vertexDegree > 0),histBins);
            
            % ----- Branch Direction Analysis ----- %
            
            [branchDir{iMov}{iFrame},branchAng{iMov}{iFrame}] = calcBranchTipDirections(currSkel.edges,currSkel.edgePaths,size(currSkel.vertices,1));                                                                                                
            
            
            % ----- Tip Path Length Analysis ---- %
                        
            [tmp1,tmp2] = analyzeSkeletonTipPaths(currSkel.vertices,currSkel.edges,currSkel.edgePaths,currSkel.edgeLabels);
            tipPathLen{iMov}{iFrame} = tmp1(~isnan(tmp1) & ~isinf(tmp1)) .* MA(iMov).pixelSize_ / 1e3;
            tipPathN{iMov}{iFrame} = cellfun(@numel,tmp2(~isnan(tmp1) & ~isinf(tmp1)));                        
            
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
        print(pOpt{:},[currOutDir filesep 'Branch Degree Distribution.eps']);
        
        
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
        
        
        % ----- Tip Path-Length Figure ---- %
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
        
    
    else
        warning('MIGRATION3D:ppSkelGraph:noPruning',...
            ['Movie ' num2str(iMov) ' does not have valid skeleton pruning - not analyzing!']);
    end
    
    if ~p.BatchMode        
        waitbar(iMov/nMovies,wtBar)
    end
        
end

close(wtBar);

% ------ Combined (All-Movie) Branch-Point-Degree Data ------- %

combMean = nanmean(meanHist,1);
combSTD = nanstd(meanHist,1);

if p.BatchMode
    cmFig(iMov) = figure('Visible','off');
else
    cmFig(iMov) = figure;
end
subplot(2,1,1)
bar(combMean)
hold on
errorbar(combMean,combSTD,'.r')
set(gca,'XTickLabel',histBins)
xlabel('Branch Point Degree')
ylabel('Average Count Per Frame')
title('Branch point degree distribution, average per-frame')
subplot(2,1,2)
combProbDist = nanmean(meanHist ./ repmat(sum(meanHist,2),[1 numel(histBins)]),1);
combProbSTD = nanstd(meanHist ./ repmat(sum(meanHist,2),[1 numel(histBins)]),[],1);
bar(combProbDist)
hold on
errorbar(combProbDist,combProbSTD,'.r')
set(gca,'XTickLabel',histBins)
xlabel('Branch Point Degree')
ylabel('Average Probability')
title('Branch point degree probability distribution, all movies')

print(pOpt{:},[p.OutputDirectory filesep 'Combined Branch Degree Distribution.eps']);


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
title({'Unaligned, combined Phi Histogram, all Movies',['Combined Phi STD: ' num2str(stdPhi)]})

combThetaHist = histc(combAngData(:,1),thetaBins);
meanTheta = nanmean(combAngData(:,1));
stdTheta = nanstd(combAngData(:,1));
subplot(2,1,2);
bar(thetaBins,combThetaHist)
xlim([-pi pi])
xlabel('Theta, Radians')
ylabel('Count')
title({'Unaligned, combined Theta Histogram, all Movies',['Combined Theta STD: ' num2str(stdTheta)]})
print(pOpt{:},[p.OutputDirectory filesep 'Combined Branch Angle Distributions.eps']);




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
%errorbar(combPathLenBins,mean(combPathHist,1),std(combPathHist,[],1),'.r')
title('Combined Tip-To-Body Path-Length Probability Distribution, All Movies')
xlabel('Path Length, microns')
ylabel('Average Probability')

subplot(2,1,2)
bar(combPathNBins,mean(combPathNHist,1))
hold on
%errorbar(combPathNBins,mean(combPathNHist,1),std(combPathNHist,[],1),'.r')
title('Combined Tip-To-Body Path Complexity Probability Distribution, All Movies')
xlabel('Path Complexity, # of Vertices')
ylabel('Average Probability')







%TEMP - in batch mode, go through and close all the figures when you're
%done!!!!

