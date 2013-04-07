function postProcess3DMovieArrayMaskTracking(MA,varargin)


%TEMP - Yet again, writing a function sloppily to get a result at the last
%minute!! Yay time management!!


%% ------------- Parameters -------------- %%

perMovDirName = 'mask__object_tracking_post_processing';%Directory name for saving results from individual movies
perMovFName = 'centroid_tracking_post_processing';
%Print options for saving figures to disk
pOpt = {'-r300',...% dpi = 300
        '-depsc2'};% use eps format
pOptTIFF = {'-r100','-dtiff'};%150 dpi for TIF format since this is usally just used for emailing to Bob!


iProcChan = 1;%So what, your mother hard-codes parameters too!!

        
%% ------------- Input ------------ %%


%Parse all the inputs
ip = inputParser;
ip.FunctionName = mfilename;
ip.addRequired('MA',@(x)(isa(x,'MovieData3D')));
ip.addParamValue('BatchMode',false,(@(x)(numel(x)==1)));
ip.addParamValue('UseTimeInterval',[],(@(x)(numel(x)==1 && x > 0)));
ip.addParamValue('OutputDirectory',[],@ischar);%Directory will be created if it doesn't exist
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

mkClrDir(p.OutputDirectory);

objDispPerFrameOutliers = cell(nMovies,1);
pixXY = nan(nMovies,1);
timeInt = nan(nMovies,1);
hasET = false(nMovies,1);
hasPS = false(nMovies,1);
hasMG = false(nMovies,1);

objTracks = cell(nMovies,1);
objDispPerFrame = cell(nMovies,1);
objDispStartEnd = nan(nMovies,1);
objDispPerFrameSum = nan(nMovies,1);
objVelPerFrame = cell(nMovies,1);
objVelPerFrameMean = nan(nMovies,1);
objVelStartEnd = nan(nMovies,1);

if p.BatchMode
    figArgs = {'Visible','off'};
else
    figArgs = {};
end

nAngBins = 10;%Number of bins for angle histograms
degHistBins = [1 3 4 5 6 7];%Skip two - think about it.
nDegHistBins = numel(degHistBins);
phiBins = linspace(-pi/2,pi/2,nAngBins);
thetaBins = linspace(-pi,pi,nAngBins);
tipPathBinSz = 5;

nTips = cell(nMovies,1);
iTips = cell(nMovies,1);

branchRad = cell(nMovies,1);

maxLag = 6;%TEMP Stupid hard-coded max lag becuase I need to get this done fast


%% --------------- Per Movie Processing --------------- %%


for iMov = 1:nMovies
    
    
    timeInt(iMov) = MA(iMov).timeInterval_;
    pixXY(iMov) = MA(iMov).pixelSize_;%We only need the XY since the object tracks have been scaled to symmetric-voxel space

    %Check that the processing was completed successfully
    iMTProc = MA(iMov).getProcessIndex('MaskObjectTrackingProcess',1,~p.BatchMode);
    
    
    
    if ~isempty(iMTProc) && MA(iMov).processes_{iMTProc}.checkChannelOutput(iProcChan);
        
        if isempty(MA(iMov).eventTimes_)        
            nFramesPerMov(iMov) = MA(iMov).nFrames_;            
        else
            nFramesPerMov(iMov) = MA(iMov).eventTimes_(2);
            hasET(iMov) = true;
        end                
        
        %Load the object tracks and convert them to physical units so we can compare across movies
        objTracks{iMov} = squeeze(MA(iMov).processes_{iMTProc}.loadChannelOutput(iProcChan)) .* pixXY(iMov);
        
        %Remove frames with bad segmentation or other problems
        objTracks{iMov} = objTracks{iMov}(1:nFramesPerMov(iMov),:);

        
        % ----- Calc Per-Movie Motion Statistics ---- %
        
        
        %Calculate the frame-to-frame displacements
        objDispPerFrame{iMov} = sqrt(sum(diff(objTracks{iMov},1,1) .^2,2));
        %And total distance traveled
        objDispPerFrameSum(iMov) = sum(objDispPerFrame{iMov});
        %And total start-to-end striaght-line distance traveled
        objDispStartEnd(iMov) = sqrt(sum(diff(objTracks{iMov}([1 end],:),1,1) .^2,2));
        %And frame-to-frame velocities in nm/s
        objVelPerFrame{iMov} = objDispPerFrame{iMov} / timeInt(iMov);
        objVelPerFrameMean(iMov) = mean(objVelPerFrame{iMov});
        %And the effective velocity from start-to-end
        objVelStartEnd(iMov) = objDispStartEnd(iMov) / (timeInt(iMov)*nFramesPerMov(iMov));
        
        
        %Check for outliers - the centermost-point method is not completely
        %independent of cell shape so can occasionally make large jumps.
        objDispPerFrameOutliers{iMov} = detectOutliers(objDispPerFrame{iMov},3);
        
        %Set-up output for this movie        
        currOutDir = [MA(iMov).outputDirectory_ filesep perMovDirName];
        mkClrDir(currOutDir);        
        
        % ----- Per- Movie Figures ----- %
        
        dispFig = figure(figArgs{:});
        dispEdges = [(1:nFramesPerMov(iMov)-1)' (2:nFramesPerMov(iMov))'];
        plot3(objTracks{iMov}(:,1),objTracks{iMov}(:,2),objTracks{iMov}(:,3),'-','MarkerSize',5);        
        hold on,axis equal
        for j = objDispPerFrameOutliers{iMov}'
            plot3(objTracks{iMov}(dispEdges(j,:),1),objTracks{iMov}(dispEdges(j,:),2),objTracks{iMov}(dispEdges(j,:),3),'r','LineWidth',2);        
        end
        if ~isempty(objDispPerFrameOutliers{iMov})
            legend('Centroid Track','Outlier Displacements')        
        end
        %Color-code the centroid locations to indicate time
        timeCols = jet(nFramesPerMov(iMov));
        scatter3(objTracks{iMov}(:,1),objTracks{iMov}(:,2),objTracks{iMov}(:,3),ones(nFramesPerMov(iMov),1)*pixXY(iMov),timeCols,'filled')
        
        xlabel('X position , nm')
        ylabel('Y position , nm')
        zlabel('Z position , nm')                                        
        title('Centroid position over time, blue is first frame red is last')       
        
        figName = [currOutDir filesep 'centroid tracks'];
        print(pOpt{:},[figName '.eps']);
        print(pOptTIFF{:},[figName '.tif']);
        hgsave([figName '.fig'])    
        
        statFig = figure(figArgs{:});
        subplot(2,1,1);
        hist(objVelPerFrame{iMov})
        title({'Per-Frame Velocity Histogram',...
              ['Total Displacement =  ' num2str(objDispStartEnd(iMov)) 'nm'],...
              ['Total Displacement / time =  ' num2str(objVelStartEnd(iMov)) 'nm/s'],...
              ['Mean F-t-F Vel =  ' num2str(objVelPerFrameMean(iMov)) 'nm/s'],...
              ['Total F-t-F disp = ' num2str(objDispPerFrameSum(iMov)) 'nm']});
        subplot(2,1,2);
        plot(0:timeInt(iMov):(nFramesPerMov(iMov)-2)*timeInt(iMov),objVelPerFrame{iMov},'.-')
        xlabel('Time, s')
        ylabel('Velocity, nm/s')
        title('Centroid Velocity vs. Time')
        figName = [currOutDir filesep 'velocity statistics'];
        print(pOpt{:},[figName '.eps']);
        print(pOptTIFF{:},[figName '.tif']);
        hgsave([figName '.fig'])    
                      
    else
        warning('MIGRATION3D:ppMaskTrack:noTracking',...
            ['Movie ' num2str(iMov) ' does not have valid mask object tracking - not analyzing!']);
    end
    if p.BatchMode
        close all
    end
    
    iSkProc = MA(iMov).getProcessIndex('SkeletonPruningProcess',1,~p.BatchMode);
    
    %TEMP - this duplicates the shit in postPRocessSKeletons, so we should
    %just save it to disk and load it here...??? Or should we just
    %cut-and-paste and fucking graduate!?
    
    if ~isempty(iSkProc) && MA(iMov).processes_{iSkProc}.checkChannelOutput(iProcChan);
        
        hasPS(iMov) = true;
        nTips{iMov} = nan(nFramesPerMov(iMov),1);
        iTips{iMov} = cell(nFramesPerMov(iMov),1);
        branchDir{iMov} = cell(nFramesPerMov(iMov),1);
        branchAng{iMov} = cell(nFramesPerMov(iMov),1);
        tipPathLen{iMov} = cell(nFramesPerMov(iMov),1);
        tipPathN{iMov} = cell(nFramesPerMov(iMov),1);
        branchRad{iMov} = cell(nFramesPerMov(iMov),1);
        
        iMgProc = MA(iMov).getProcessIndex('MaskGeometry3DProcess',1,~p.BatchMode);
        
        if ~isempty(iMgProc) && MA(iMov).processes_{iMgProc}.checkChannelOutput(iProcChan)
            hasMG(iMov) = true;
        end
        
        for iFrame = 1:nFramesPerMov(iMov);
        
            currSkel = MA(iMov).processes_{iSkProc}.loadChannelOutput(iProcChan,iFrame);
            
            if hasMG(iMov)
                currMaskProp = MA(iMov).processes_{iMgProc}.loadChannelOutput(iProcChan,iFrame);
                branchRad{iMov}{iFrame} = branchRadii(currSkel,currMaskProp);
                branchRad{iMov}{iFrame} = cellfun(@(x)(x .* pixXY(iMov)),branchRad{iMov}{iFrame},'Unif',0);%Convert to nm
            end
            
            % ----- Vertex Degree Analysis ----- %
            
            [iTips{iMov}{iFrame},iTipEdges{iMov}{iFrame}] = findTips(currSkel.edges,size(currSkel.vertices,1));
            nTips{iMov}(iFrame) = numel(iTips{iMov}{iFrame});
            
            skelStats.vertexDegree = graphVertDegree(currSkel.edges,size(currSkel.vertices,1));             
            vdHists{iMov}(iFrame,:) = histc(skelStats.vertexDegree(skelStats.vertexDegree > 0),degHistBins);
            
            isBranch = currSkel.edgeLabels == 1;
            
            meanRadPerSkelementPerFrame{iMov}{iFrame} = cellfun(@mean,branchRad{iMov}{iFrame});
            
            meanRadPerBranchPerFrame{iMov}{iFrame} = cellfun(@mean,branchRad{iMov}{iFrame}(isBranch));
            meanRadPerFrameBranchWeighted{iMov}(iFrame) = mean(meanRadPerBranchPerFrame{iMov}{iFrame});%Averages the radius with each branch contributing equally
            meanRadPerFramePointWeighted{iMov}(iFrame) = mean(vertcat(branchRad{iMov}{iFrame}{isBranch}));%Averages the radius with each POINT contributing equally, so that longer branches contribute more
            meanTipRadPerBranchPerFrame{iMov}{iFrame} = cellfun(@mean,branchRad{iMov}{iFrame}(iTipEdges{iMov}{iFrame}));
            meanTipRadPerFrameBranchWeighted{iMov}(iFrame) = mean(meanRadPerSkelementPerFrame{iMov}{iFrame}(iTipEdges{iMov}{iFrame}));%Averages the radius with each branch contributing equally
            meanTipRadPerFramePointWeighted{iMov}(iFrame) = mean(vertcat(branchRad{iMov}{iFrame}{(iTipEdges{iMov}{iFrame})}));%Averages the radius with each POINT contributing equally, so that longer branches contribute more
                        
            totalMeanRadPerFrameBranchWeighted{iMov}(iFrame) = sum(meanRadPerBranchPerFrame{iMov}{iFrame});%Sums the mean radii with each branch contributing equally
            totalMeanRadPerFramePointWeighted{iMov}(iFrame) = sum(vertcat(branchRad{iMov}{iFrame}{isBranch}));%Sums the radius with each POINT contributing equally, so that longer branches contribute more
            %totalMeanTipRadPerBranchPerFrame{iMov}{iFrame} = cellfun(@sum,branchRad{iMov}{iFrame}(iTipEdges{iMov}{iFrame}));
            totalMeanTipRadPerFrameBranchWeighted{iMov}(iFrame) = sum(meanRadPerSkelementPerFrame{iMov}{iFrame}(iTipEdges{iMov}{iFrame}));%Averages the radius with each branch contributing equally
            totalMeanTipRadPerFramePointWeighted{iMov}(iFrame) = sum(vertcat(branchRad{iMov}{iFrame}{(iTipEdges{iMov}{iFrame})}));%Averages the radius with each POINT contributing equally, so that longer branches contribute more
            
            % ----- Branch Direction Analysis ----- %
            
            [branchDir{iMov}{iFrame},branchAng{iMov}{iFrame}] = calcBranchTipDirections(currSkel.edges,currSkel.edgePaths,size(currSkel.vertices,1));
            
            
            % ----- Tip Path Length Analysis ---- %
                        
            [tmp1,tmp2] = analyzeSkeletonTipPaths(currSkel.vertices,currSkel.edges,currSkel.edgePaths,currSkel.edgeLabels);
            tipPathLen{iMov}{iFrame} = tmp1(~isnan(tmp1) & ~isinf(tmp1)) .* MA(iMov).pixelSize_ / 1e3;
            tipPathN{iMov}{iFrame} = cellfun(@numel,tmp2(~isnan(tmp1) & ~isinf(tmp1)));                        
            
        end        
        
        meanRadPerMovBranchWeighted(iMov) = mean(meanRadPerFrameBranchWeighted{iMov});
        meanRadPerMovPointWeighted(iMov) = mean(meanRadPerFramePointWeighted{iMov});
        meanTipRadPerMovBranchWeighted(iMov) = mean(meanTipRadPerFrameBranchWeighted{iMov});
        meanTipRadPerMovPointWeighted(iMov) = mean(meanTipRadPerFramePointWeighted{iMov});
        
        meanTotalRadPerMovBranchWeighted(iMov) = mean(totalMeanRadPerFrameBranchWeighted{iMov});
        meanTotalRadPerMovPointWeighted(iMov) = mean(totalMeanRadPerFramePointWeighted{iMov});
        meanTotalTipRadPerMovBranchWeighted(iMov) = mean(totalMeanTipRadPerFrameBranchWeighted{iMov});
        meanTotalTipRadPerMovPointWeighted(iMov) = mean(totalMeanTipRadPerFramePointWeighted{iMov});
        
        
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
        tipNMean(iMov) = mean(allTipPathN{iMov});
        tipNMed(iMov) = median(allTipPathN{iMov});
        tipLenMean(iMov) = mean(allTipPathLen{iMov});
        tipLenMed(iMov) = median(allTipPathLen{iMov});
        
        if nFramesPerMov(iMov) > 1
        
            vandbFig = figure(figArgs{:});                
            subplot(2,1,1);
            title('Per-Frame Velocity and Branch Number')
            %We have displacements for n-1 frames, so we arbitratily leave of
            %the last frame of branch number
            [axHan,hY1,hY2] = plotyy(0:timeInt(iMov):(nFramesPerMov(iMov)-2)*timeInt(iMov),objVelPerFrame{iMov}',...
                    0:timeInt(iMov):(nFramesPerMov(iMov)-2)*timeInt(iMov),nTips{iMov}(1:end-1));
            xlabel('Time, s')    
            set(get(axHan(1),'Ylabel'),'String','Instantaneous Velocity, nm/s') 
            set(get(axHan(2),'Ylabel'),'String','Tip #') 
            subplot(2,1,2);                
            [axHan,hY1,hY2] = plotyy(0:timeInt(iMov):(nFramesPerMov(iMov)-2)*timeInt(iMov),objVelPerFrame{iMov}',...
                    0:timeInt(iMov):(nFramesPerMov(iMov)-2)*timeInt(iMov),totalMeanTipRadPerFrameBranchWeighted{iMov}(1:end-1));
            xlabel('Time, s')    
            set(get(axHan(1),'Ylabel'),'String','Instantaneous Velocity, nm/s') 
            set(get(axHan(2),'Ylabel'),'String','Radius Weighted Tip #') 

            figName = [currOutDir filesep 'velocity and branch statistics time series'];
            print(pOpt{:},[figName '.eps']);
            print(pOptTIFF{:},[figName '.tif']);
            hgsave([figName '.fig'])    
        end
        
        % ---- Per-Movie Cross-Corr Figures ---- %%
        
        if nFramesPerMov(iMov) >= maxLag+1
            tLags = -(maxLag*timeInt(iMov)):timeInt(iMov):(maxLag*timeInt(iMov));

            [tmp,~,tmpB] = crosscorr(totalMeanTipRadPerFrameBranchWeighted{iMov}(1:end-1)',objVelPerFrame{iMov},maxLag);
            figure;
            plot(tLags,tmp)
            hold on
            plot(xlim,ones(1,2)*tmpB(1),'--r')
            plot(xlim,ones(1,2)*tmpB(2),'--r')
            plot([0 0 ],ylim,'--k')
            plot(xlim,[0 0 ],'--k')
            xlabel('Delay, seconds (Positive Meanse Velocity follows branching')
            ylabel('Cross Correlation')
            title('Cross-Corr, Radius-Weighted Branch Number and Instantaneous Velocity')

            figName = [currOutDir filesep 'cross corr velocity and radius weighted branch number'];
            print(pOpt{:},[figName '.eps']);
            print(pOptTIFF{:},[figName '.tif']);
            hgsave([figName '.fig'])    

            ccPerCellRadWtVel(iMov,:) = tmp;
            cbPerCellRadWtVel(iMov,:) = tmpB;
        
        else
            ccPerCellRadWtVel(iMov,:) = NaN;
            cbPerCellRadWtVel(iMov,:) = NaN;
        end
        if nFramesPerMov(iMov) >= maxLag+1

            [tmp,~,tmpB] = crosscorr(nTips{iMov}(1:end-1),objVelPerFrame{iMov},maxLag);
            figure;
            plot(tLags,tmp(:,1))
            hold on
            plot(xlim,ones(1,2)*tmpB(1),'--r')
            plot(xlim,ones(1,2)*tmpB(2),'--r')
            plot([0 0 ],ylim,'--k')
            plot(xlim,[0 0 ],'--k')
            xlabel('Delay, seconds (Positive means velocity follows branching')
            ylabel('Cross Correlation')
            title('Cross-Corr, Thresholded Tip-Count and Instantaneous Velocity')

            figName = [currOutDir filesep 'cross corr velocity and thresholded tip count'];
            print(pOpt{:},[figName '.eps']);
            print(pOptTIFF{:},[figName '.tif']);
            hgsave([figName '.fig'])    

            ccPerCellThreshTipVel(iMov,:) = tmp;
            cbPerCellThreshTipVel(iMov,:) = tmpB;
        else            
            ccPerCellThreshTipVel(iMov,:) = NaN;
            cbPerCellThreshTipVel(iMov,:) = NaN;            
        end
        
%         [tmp,tmpA] = modifiedKendallCorr(nTips{iMov}(1:end-1),objVelPerFrame{iMov},maxLag-2,.05,true,maxLag);
%         figure;
%         plot(tLags,tmp)
%         hold on
%         %plot(xlim,1.96/sqrt(nFramesPerMov(iMov))*ones(1,2),'--r')
%         %plot(xlim,-1.96/sqrt(nFramesPerMov(iMov))*ones(1,2),'--r')
%         plot([0 0 ],ylim,'--k')
%         plot(xlim,[0 0 ],'--k')
%         xlabel('Delay, seconds (Positive means velocity follows branching')
%         ylabel('Cross Correlation')
%         title('Modified Kendall Corr, Thresholded Tip-Count and Instantaneous Velocity')
%         
%         kcPerCellThreshTipVel(iMov,:) = tmp;
%         kcPerCellThreshTipVelA(iMov,:) = tmpA;
%         figName = [currOutDir filesep 'kendall corr corr velocity and thresholded tip count'];
%         print(pOpt{:},[figName '.eps']);
%         print(pOptTIFF{:},[figName '.tif']);
%         hgsave([figName '.fig'])    

        
            
    end

end

disp([num2str(nnz(hasET)) ' of ' num2str(nMovies) ' had event times and were time-truncated.'])
if nnz(hasPS) ~= nMovies
    warning('Not all movies had skeleton post-processing!!!')
end
if nnz(hasMG) ~= nMovies
    warning('Not all movies had mask geometry analysis!!')
end

%% ------------- All-Movie Combined Analysis ----------- %%

velBins = linspace(0,150,20);
allVelPerFrame = vertcat(objVelPerFrame{:});

% ---- Velocity Histograms ------ %


vHistFitPF = figure(figArgs{:});
subplot(1,4,1)
hist(allVelPerFrame,velBins)
title({'Instantaneous Velocity Distribution',...
        ['All Cells, n=' num2str(nMovies)],...
        ['Combined Mean = ' num2str(mean(allVelPerFrame)) 'nm/s'],...
        ['Robust Mean = ' num2str(robustMean(allVelPerFrame)) 'nm/s']});
xlabel('nm/s')
ylabel('Total Frame Count')

subplot(1,4,2)
hist(objDispPerFrameSum)
title({'Total Path-Length',...
        ['All Cells, n=' num2str(nMovies)],...
        ['Combined Mean = ' num2str(mean(objDispPerFrameSum)) 'nm'],...
        ['Robust Mean = ' num2str(robustMean(objDispPerFrameSum)) 'nm']});
xlabel('nm')
ylabel('# Cells')

subplot(1,4,3)
hist(objVelStartEnd)
title({'Start-End Velocity',...
        ['All Cells, n=' num2str(nMovies)],...
        ['Combined Mean = ' num2str(mean(objVelStartEnd)) 'nm/s'],...
        ['Robust Mean = ' num2str(robustMean(objVelStartEnd)) 'nm/s']});
xlabel('nm/s')
ylabel('# Cells')

subplot(1,4,4)
hist(objDispStartEnd)
title({'Start-End Displacement',...
        ['All Cells, n=' num2str(nMovies)],...
        ['Combined Mean = ' num2str(mean(objDispStartEnd)) 'nm'],...
        ['Robust Mean = ' num2str(robustMean(objDispStartEnd)) 'nm']});
xlabel('nm')
ylabel('# Cells')


figName = [p.OutputDirectory filesep 'velocity histograms'];
print(pOpt{:},[figName '.eps']);
print(pOptTIFF{:},[figName '.tif']);
hgsave([figName '.fig'])    



%% ------------ Velocity / Branch Structure All-Movie Correlation ------ %%


% ------ Instantaneous Velocity vs. X Figures - Single Corr ALl Cells ----- %

if ~isempty(p.UseTimeInterval)
    useTime  = p.UseTimeInterval;
else
    %Otherwise just use the most common time interval
    useTime = mode(timeInt);
end

isCorrTimeInt = timeInt == useTime;

padVelStruc(1:nMovies,1) = struct('observations',[]);

for j = 1:nMovies
    %We don't have velocities for the last frame, so we pad with NaN to make the
    %plotting/correlation easier
    padVelPerFrame{j} = vertcat(objVelPerFrame{j},NaN);    
    if isCorrTimeInt(j)
        %Put them in a struc for cross-corr analysis
        padVelStruc(j).observations = padVelPerFrame{j};
        nTipsStruc(j).observations = nTips{j};
        totMeanPointWeightRadStruc(j).observations = totalMeanRadPerFramePointWeighted{j}';
        totMeanBranchWeightRadStruc(j).observations = totalMeanRadPerFrameBranchWeighted{j}';
        totMeanTipPointWeightRadStruc(j).observations = totalMeanTipRadPerFramePointWeighted{j}';
        totMeanTipBranchWeightRadStruc(j).observations = totalMeanTipRadPerFrameBranchWeighted{j}';
    end
end

padVelStruc = padVelStruc(isCorrTimeInt);
nTipsStruc = nTipsStruc(isCorrTimeInt);
totMeanPointWeightRadStruc = totMeanPointWeightRadStruc(isCorrTimeInt);
totMeanBranchWeightRadStruc = totMeanBranchWeightRadStruc(isCorrTimeInt);
totMeanTipPointWeightRadStruc = totMeanTipPointWeightRadStruc(isCorrTimeInt);
totMeanTipBranchWeightRadStruc = totMeanTipBranchWeightRadStruc(isCorrTimeInt);


velTipXcorr = crossCorr(padVelStruc,nTipsStruc,maxLag);

allPadVelPerFrame = vertcat(padVelPerFrame{:});
allNTipsPerFrame = vertcat(nTips{:});

justTheTipFig = figure(figArgs{:});
subplot(1,2,1);
plot(allPadVelPerFrame,allNTipsPerFrame,'.')
title('Scatter, Instantaneous velocity and Tip Count')
xlabel('Instantaneous Velocity, nm/s')
ylabel('Number of Branch Tips')

subplot(1,2,2);
plot(-useTime*maxLag:useTime:useTime*maxLag,velTipXcorr(:,1),'.-')
hold on
plot(xlim,ones(1,2)*1.96/sqrt(sum(nFramesPerMov(isCorrTimeInt))),'--')
plot(xlim,ones(1,2)*-1.96/sqrt(sum(nFramesPerMov(isCorrTimeInt))),'--')
plot(xlim,[0 0],'--k')
plot([0 0],ylim,'--k')
title({'Cross Corr, Tip Count and Instantaneous Velocity',...
        'Positive Delay means nTips Follows Velocity',...
        [num2str(useTime) 's data only']})
xlabel('Time Delay, s')
ylabel('Correlation')

figName = [p.OutputDirectory filesep 'instantaneous velocity and tip number'];
print(pOpt{:},[figName '.eps']);
print(pOptTIFF{:},[figName '.tif']);
hgsave([figName '.fig'])    


radFig = figure(figArgs{:});
subplot(1,2,1)
plot(allPadVelPerFrame,horzcat(totalMeanRadPerFramePointWeighted{:}),'.')
xlabel('Instantaneous Velocity')
ylabel('Point-Weighted Branch Radius (~branch volume) per frame, nm')
title('Point-Weighted Total Branch Radius (~Branch Volume) vs. Instantaneous Velocity')     
subplot(1,2,2)
velTipXcorr = crossCorr(padVelStruc,totMeanPointWeightRadStruc,maxLag);
plot(-useTime*maxLag:useTime:useTime*maxLag,velTipXcorr(:,1),'.-')
hold on
plot(xlim,ones(1,2)*-1.96/sqrt(sum(nFramesPerMov(isCorrTimeInt))),'--')
plot(xlim,ones(1,2)*1.96/sqrt(sum(nFramesPerMov(isCorrTimeInt))),'--')
xlabel('time delay, s')
ylabel('correlation')
  
figName = [p.OutputDirectory filesep 'instantaneous velocity and point weighted mean branch radius'];
print(pOpt{:},[figName '.eps']);
print(pOptTIFF{:},[figName '.tif']);
hgsave([figName '.fig'])    


radFig = figure(figArgs{:});
subplot(1,2,1)
plot(allPadVelPerFrame,horzcat(totalMeanRadPerFrameBranchWeighted{:}),'.')
xlabel('Instantaneous Velocity')
ylabel('Radius-weighted skeleton element number per frame, nm')
title('Radius-weighted skeleton element number vs. Instantaneous Velocity')     
subplot(1,2,2)
velTipXcorr = crossCorr(padVelStruc,totMeanBranchWeightRadStruc,maxLag);
plot(-useTime*maxLag:useTime:useTime*maxLag,velTipXcorr(:,1),'.-')
hold on
plot(xlim,ones(1,2)*-1.96/sqrt(sum(nFramesPerMov(isCorrTimeInt))),'--')
plot(xlim,ones(1,2)*1.96/sqrt(sum(nFramesPerMov(isCorrTimeInt))),'--')
xlabel('time delay, s')
ylabel('correlation')

figName = [p.OutputDirectory filesep 'instantaneous velocity and radius weighted mean branch number'];
print(pOpt{:},[figName '.eps']);
print(pOptTIFF{:},[figName '.tif']);
hgsave([figName '.fig'])    


radFig = figure(figArgs{:});
subplot(1,2,1)
plot(allPadVelPerFrame,horzcat(totalMeanTipRadPerFrameBranchWeighted{:}),'.')
xlabel('Instantaneous Velocity')
ylabel('Radius-Weighted Tip Number per frame, nm')
title('Radius-Weighted Tip Number vs. Instantaneous Velocity')     
subplot(1,2,2)
velTipXcorr = crossCorr(padVelStruc,totMeanTipBranchWeightRadStruc,maxLag);
plot(-useTime*maxLag:useTime:useTime*maxLag,velTipXcorr(:,1),'.-')
hold on
plot(-useTime*maxLag:useTime:useTime*maxLag,velTipXcorr(:,1) + velTipXcorr(:,2),'--')
plot(-useTime*maxLag:useTime:useTime*maxLag,velTipXcorr(:,1) - velTipXcorr(:,2),'--')
plot(xlim,[0 0],'--k')
plot([0 0],ylim,'--k')
% plot(xlim,ones(1,2)*-1.96/sqrt(sum(nFramesPerMov(isCorrTimeInt))),'--')
% plot(xlim,ones(1,2)*1.96/sqrt(sum(nFramesPerMov(isCorrTimeInt))),'--')
xlabel('time delay, s')
ylabel('correlation')
title({'Cross Corr, Radius-Weighted Tip Count and Instantaneous Velocity',...
        'Positive Delay means tip # Follows Velocity',...
        ['n=' num2str(nnz(isCorrTimeInt)) ' cells, ' num2str(useTime) 's data only']})
figName = [p.OutputDirectory filesep 'instantaneous velocity and radius weighted mean tip number'];
print(pOpt{:},[figName '.eps']);
print(pOptTIFF{:},[figName '.tif']);
hgsave([figName '.fig'])    


radFig = figure(figArgs{:});
subplot(1,2,1)
plot(allPadVelPerFrame,horzcat(totalMeanTipRadPerFramePointWeighted{:}),'.')
xlabel('Instantaneous Velocity')
ylabel('Point-Weighted Total Tip Radius per frame, nm')
title('Point-Weighted Total Tip Radius vs. Instantaneous Velocity')     
subplot(1,2,2)
velTipXcorr = crossCorr(padVelStruc,totMeanTipPointWeightRadStruc,maxLag);
plot(-useTime*maxLag:useTime:useTime*maxLag,velTipXcorr(:,1),'.-')
hold on
plot(xlim,ones(1,2)*-1.96/sqrt(sum(nFramesPerMov(isCorrTimeInt))),'--')
plot(xlim,ones(1,2)*1.96/sqrt(sum(nFramesPerMov(isCorrTimeInt))),'--')
xlabel('time delay, s')
ylabel('correlation')

figName = [p.OutputDirectory filesep 'instantaneous velocity and point weighted mean tip radius'];
print(pOpt{:},[figName '.eps']);
print(pOptTIFF{:},[figName '.tif']);
hgsave([figName '.fig'])    

% ------ Instantaneous Velocity vs. X Figures - Average of Per-Cell Corrs ----- %



[meanOfPerCellCCradWt,ciOfPerCellCCradWt] = correlationBootstrap(ccPerCellRadWtVel',cbPerCellRadWtVel(:,1)');

figure
plot(-useTime*maxLag:useTime:useTime*maxLag,meanOfPerCellCCradWt)
hold on
plot(-useTime*maxLag:useTime:useTime*maxLag,ciOfPerCellCCradWt(1,:),'--r')
legend('Correlation','Boostrapped 95% C.I.')
plot(-useTime*maxLag:useTime:useTime*maxLag,ciOfPerCellCCradWt(2,:),'--r')
plot(xlim,[ 0 0 ],'--k')
plot([ 0 0 ],ylim,'--k')
xlabel('Time Delay, s')
ylabel('Correlation')
title({'Radius Weighted Tip Count and Instantaneous Velocity',...
    ['Average of Per-Cell Correlations, n=' num2str(nMovies)]})

figName = [p.OutputDirectory filesep 'instantaneous velocity and weighted mean tip number avg of per-cell corr'];
print(pOpt{:},[figName '.eps']);
print(pOptTIFF{:},[figName '.tif']);
hgsave([figName '.fig'])    

[meanOfPerCellCCThreshTip,ciOfPerCellCCThreshTip] = correlationBootstrap(ccPerCellThreshTipVel',cbPerCellThreshTipVel(:,1)');

figure
plot(-useTime*maxLag:useTime:useTime*maxLag,meanOfPerCellCCThreshTip)
hold on
plot(-useTime*maxLag:useTime:useTime*maxLag,ciOfPerCellCCThreshTip(1,:),'--r')
legend('Correlation','Boostrapped 95% C.I.')
plot(-useTime*maxLag:useTime:useTime*maxLag,ciOfPerCellCCThreshTip(2,:),'--r')
plot(xlim,[ 0 0 ],'--k')
plot([ 0 0 ],ylim,'--k')
xlabel('Time Delay, s')
ylabel('Correlation')
title({'Thresholded Tip Count and Instantaneous Velocity',...
    ['Average of Per-Cell Correlations, n=' num2str(nMovies)]})

figName = [p.OutputDirectory filesep 'instantaneous velocity and thresholded tip number avg of per-cell corr'];
print(pOpt{:},[figName '.eps']);
print(pOptTIFF{:},[figName '.tif']);
hgsave([figName '.fig'])    



%% ------------ Mean  Velocity vs. X Figures ------------ %

meanTipsPerMov = cellfun(@mean,nTips);
stdTipsPerMov = cellfun(@std,nTips);

justTheTipFig = figure(figArgs{:});
subplot(1,2,1);
x = objVelPerFrameMean;
y = meanTipsPerMov;
hasOb = ~(isnan(x) | isnan(y));
x = x(hasOb);
y = y(hasOb);
plot(x,y,'.')
hold on
[fitObj,gofStats] = fit(x,y,'poly1');
yFit = feval(fitObj,x);
fitCI = confint(fitObj);
plot(x,yFit,'r')
xlabel('Mean Per-Frame Velocity, nm/s')
ylabel('Average # Tips Per Frame')
legend('Value','Fit')
title(['m=' num2str(fitObj.p1) ', 95%CI= ' num2str(fitCI(1,1)) ' to ' num2str(fitCI(2,1))])
subplot(1,2,2);
x = objVelStartEnd;
y = meanTipsPerMov;
plot(x,y,'.')
hold on
[fitObj,gofStats] = fit(x,y,'poly1');
yFit = feval(fitObj,x);
fitCI = confint(fitObj);
plot(x,yFit,'r')
plot(x,y,'.')
hold on
xlabel('Straight-Line Displacement/Time, nm/s')
ylabel('Average # Tips Per Frame')
legend('Value','Fit')
title(['m=' num2str(fitObj.p1) ', 95%CI= ' num2str(fitCI(1,1)) ' to ' num2str(fitCI(2,1))])

figName = [p.OutputDirectory filesep 'mean velocity and mean tip number'];
print(pOpt{:},[figName '.eps']);
print(pOptTIFF{:},[figName '.tif']);
hgsave([figName '.fig'])    


figure(figArgs{:})
subplot(1,2,1);
x = objVelPerFrameMean;
y = meanTotalRadPerMovBranchWeighted(:);
hasOb = ~(isnan(x) | isnan(y));
x = x(hasOb);
y = y(hasOb);
plot(x,y,'.')
hold on
[fitObj,gofStats] = fit(x,y,'poly1');
yFit = feval(fitObj,x);
fitCI = confint(fitObj);
plot(x,yFit,'r')
xlabel('Mean Per-Frame Velocity, nm/s')
ylabel('Radius-Weighted Branch # Per Frame')
subplot(1,2,2);
x = objVelStartEnd;
y = meanTotalRadPerMovBranchWeighted(:);
plot(x,y,'.')
hold on
[fitObj,gofStats] = fit(x,y,'poly1');
yFit = feval(fitObj,x);
fitCI = confint(fitObj);
plot(x,yFit,'r')
title(['m=' num2str(fitObj.p1) ', 95%CI= ' num2str(fitCI(1,1)) ' to ' num2str(fitCI(2,1))])
xlabel('Straight-Line Displacement/Time, nm/s')
ylabel('Radius-Weighted Branch # Per Frame')

figName = [p.OutputDirectory filesep 'mean velocity and mean radius weighted tip number'];
print(pOpt{:},[figName '.eps']);
print(pOptTIFF{:},[figName '.tif']);
hgsave([figName '.fig'])    


figure(figArgs{:})
subplot(1,2,1);
x = objVelPerFrameMean;
y = meanTotalRadPerMovPointWeighted(:);
hasOb = ~(isnan(x) | isnan(y));
x = x(hasOb);
y = y(hasOb);
plot(x,y,'.')
hold on
[fitObj,gofStats] = fit(x,y,'poly1');
yFit = feval(fitObj,x);
fitCI = confint(fitObj);
plot(x,yFit,'r')
title(['m=' num2str(fitObj.p1) ', 95%CI= ' num2str(fitCI(1,1)) ' to ' num2str(fitCI(2,1))])
xlabel('Mean Per-Frame Velocity, nm/s')
ylabel('Point-Weighted Total Branch Radius Per Frame')
subplot(1,2,2);
x = objVelStartEnd;
y = meanTotalRadPerMovPointWeighted(:);
plot(x,y,'.')
hold on
[fitObj,gofStats] = fit(x,y,'poly1');
yFit = feval(fitObj,x);
fitCI = confint(fitObj);
plot(x,yFit,'r')
title(['m=' num2str(fitObj.p1) ', 95%CI= ' num2str(fitCI(1,1)) ' to ' num2str(fitCI(2,1))])
xlabel('Straight-Line Displacement/Time, nm/s')
ylabel('Point-Weighted Total Branch Radius Per Frame')

figName = [p.OutputDirectory filesep 'mean velocity and mean total point-weighted branch radius '];
print(pOpt{:},[figName '.eps']);
print(pOptTIFF{:},[figName '.tif']);
hgsave([figName '.fig'])    


figure(figArgs{:})
subplot(1,2,1);
x = objVelPerFrameMean;
y = meanTotalTipRadPerMovBranchWeighted(:);
hasOb = ~(isnan(x) | isnan(y));
x = x(hasOb);
y = y(hasOb);
plot(x,y,'.')
hold on
[fitObj,gofStats] = fit(x,y,'poly1');
yFit = feval(fitObj,x);
fitCI = confint(fitObj);
plot(x,yFit,'r')
title(['m=' num2str(fitObj.p1) ', 95%CI= ' num2str(fitCI(1,1)) ' to ' num2str(fitCI(2,1))])
xlabel('Mean Per-Frame Velocity, nm/s')
ylabel('Point-Weighted Total Tip Radius Per Frame')
subplot(1,2,2);
x = objVelStartEnd;
y = meanTotalTipRadPerMovBranchWeighted(:);
plot(x,y,'.')
hold on
[fitObj,gofStats] = fit(x,y,'poly1');
yFit = feval(fitObj,x);
fitCI = confint(fitObj);
plot(x,yFit,'r')
title(['m=' num2str(fitObj.p1) ', 95%CI= ' num2str(fitCI(1,1)) ' to ' num2str(fitCI(2,1))])
xlabel('Straight-Line Displacement/Time, nm/s')
ylabel('Point-Weighted Total Tip Radius Per Frame')

figName = [p.OutputDirectory filesep 'mean velocity and mean total point-weighted tip radius '];
print(pOpt{:},[figName '.eps']);
print(pOptTIFF{:},[figName '.tif']);
hgsave([figName '.fig'])    

% ------ Tip Path Vs. Velocity ---- %%

%Complexity
figure(figArgs{:})
subplot(1,2,1);
x = objVelPerFrameMean;
y = tipNMean(:);
hasOb = ~(isnan(x) | isnan(y));
x = x(hasOb);
y = y(hasOb);
plot(x,y,'.')
hold on
[fitObj,gofStats] = fit(x,y,'poly1');
yFit = feval(fitObj,x);
fitCI = confint(fitObj);
plot(x,yFit,'r')
title(['m=' num2str(fitObj.p1) ', 95%CI= ' num2str(fitCI(1,1)) ' to ' num2str(fitCI(2,1))])
xlabel('Mean Per-Frame Velocity, nm/s')
ylabel('Mean Tip Path Complexity All Frames, # of Vertices')
subplot(1,2,2);
x = objVelStartEnd;
y = tipNMean(:);
plot(x,y,'.')
hold on
[fitObj,gofStats] = fit(x,y,'poly1');
yFit = feval(fitObj,x);
fitCI = confint(fitObj);
plot(x,yFit,'r')
title(['m=' num2str(fitObj.p1) ', 95%CI= ' num2str(fitCI(1,1)) ' to ' num2str(fitCI(2,1))])
xlabel('Straight-Line Displacement/Time, nm/s')
ylabel('Mean Tip Path Complexity All Frames, # of Vertices')

figName = [p.OutputDirectory filesep 'mean velocity and mean tip path complexity'];
print(pOpt{:},[figName '.eps']);
print(pOptTIFF{:},[figName '.tif']);
hgsave([figName '.fig'])    

%Length
figure(figArgs{:})
subplot(1,2,1);
x = objVelPerFrameMean;
y = tipLenMean(:);
hasOb = ~(isnan(x) | isnan(y));
x = x(hasOb);
y = y(hasOb);
plot(x,y,'.')
hold on
[fitObj,gofStats] = fit(x,y,'poly1');
yFit = feval(fitObj,x);
fitCI = confint(fitObj);
plot(x,yFit,'r')
title(['m=' num2str(fitObj.p1) ', 95%CI= ' num2str(fitCI(1,1)) ' to ' num2str(fitCI(2,1))])
xlabel('Mean Per-Frame Velocity, nm/s')
ylabel('Mean Tip Path Length All Frames, microns')
subplot(1,2,2);
x = objVelStartEnd;
y = tipLenMean(:);
plot(x,y,'.')
hold on
[fitObj,gofStats] = fit(x,y,'poly1');
yFit = feval(fitObj,x);
fitCI = confint(fitObj);
plot(x,yFit,'r')
title(['m=' num2str(fitObj.p1) ', 95%CI= ' num2str(fitCI(1,1)) ' to ' num2str(fitCI(2,1))])
xlabel('Straight-Line Displacement/Time, nm/s')
ylabel('Mean Tip Path Length All Frames, microns')

figName = [p.OutputDirectory filesep 'mean velocity and mean tip path length'];
print(pOpt{:},[figName '.eps']);
print(pOptTIFF{:},[figName '.tif']);
hgsave([figName '.fig'])    

% ------ Velocity AutoCorr ---- %%
 
ivAc = autoCorr(padVelStruc,maxLag);
figure(figArgs{:})
plot(0:useTime:useTime*maxLag,ivAc(:,1),'.-')
hold on
title({'Autocorrelation of instantaneous velocity',...        
        ['n=' num2str(nnz(isCorrTimeInt)) ' cells, ' num2str(useTime) 's data only']})
plot(xlim,ones(1,2)*-1.96/sqrt(sum(nFramesPerMov(isCorrTimeInt))),'--')
plot(xlim,ones(1,2)*+1.96/sqrt(sum(nFramesPerMov(isCorrTimeInt))),'--')
xlabel('Time Delay, s')
ylabel('Autocorrelation')

figName = [p.OutputDirectory filesep 'instantaneous velocity autocorrelation'];
print(pOpt{:},[figName '.eps']);
print(pOptTIFF{:},[figName '.tif']);
hgsave([figName '.fig'])    




jkl=1;




