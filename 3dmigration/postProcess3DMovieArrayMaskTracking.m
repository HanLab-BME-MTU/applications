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
ip.addParamValue('BranchAngleDisplacementThresh',.5e4,(@(x)(numel(x)==1 && x > 0)));%Displacement threshold (in nm) for getting branch angles only of cells that actually  move
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
objDispVecStartEnd = nan(nMovies,3);
objDispVecStartEndUnit = nan(nMovies,3);
objDispPerFrameSum = nan(nMovies,1);
objVelPerFrame = cell(nMovies,1);
objVelPerFrameMean = nan(nMovies,1);
objVelStartEnd = nan(nMovies,1);

if p.BatchMode
    figArgs = {'Visible','off'};
else
    figArgs = {};
end

nAngBins = 20;%Number of bins for angle histograms
degHistBins = [1 3 4 5 6 7];%Skip two - think about it.
nDegHistBins = numel(degHistBins);
phiBins = linspace(-pi/2,pi/2,nAngBins);
thetaBins = linspace(-pi,pi,nAngBins);
angBins = linspace(0,pi,nAngBins);%And bins for relative angles
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
        
        tData{iMov} = 0:timeInt(iMov):(nFramesPerMov(iMov)-1)*timeInt(iMov);
        
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
        %And associated vector and unit vecor
        objDispVecStartEnd(iMov,:) = diff(objTracks{iMov}([1 end],:),1,1);
        objDispVecStartEndUnit(iMov,:) = objDispVecStartEnd(iMov,:) ./ norm(objDispVecStartEnd(iMov,:));
        [objDispVecStartEndAngle(iMov,1),objDispVecStartEndAngle(iMov,2)] = cart2sph(objDispVecStartEndUnit(iMov,1),objDispVecStartEndUnit(iMov,2),objDispVecStartEndUnit(iMov,3));
        %And frame-to-frame velocities in nm/s (yeah I know it's actually
        %speeed not velocity)
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
        branchTipDir{iMov} = cell(nFramesPerMov(iMov),1);
        branchTipAng{iMov} = cell(nFramesPerMov(iMov),1);
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
            vertDegreeMean{iMov}(iFrame) = mean(skelStats.vertexDegree);
            vertDegreeMedian{iMov}(iFrame) = median(skelStats.vertexDegree);
            vertBranchPointMean{iMov}(iFrame) = mean(skelStats.vertexDegree(skelStats.vertexDegree>2));
            vertBranchPointTotal{iMov}(iFrame) = sum(skelStats.vertexDegree(skelStats.vertexDegree>2));

                        
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
            
            % ----- Tip Path Length Analysis ---- %
                        
            [tmp1,tmp2] = analyzeSkeletonTipPaths(currSkel.vertices,currSkel.edges,currSkel.edgePaths,currSkel.edgeLabels);
            tipPathLen{iMov}{iFrame} = tmp1(~isnan(tmp1) & ~isinf(tmp1)) .* MA(iMov).pixelSize_ / 1e3;
            tipPathN{iMov}{iFrame} = cellfun(@numel,tmp2(~isnan(tmp1) & ~isinf(tmp1)));                        
            %Also get distances to centermost point 
            
            
            % ----- Branch Direction Analysis ----- %
            
            %Get Tip absolute and relative directions
            [branchTipDir{iMov}{iFrame},branchTipAng{iMov}{iFrame}] = calcBranchTipDirections(currSkel.edges,currSkel.edgePaths,size(currSkel.vertices,1));
            %And branch angles relative to net displacement vector
            branchTipAngToDisp{iMov}{iFrame} = acos(sum(branchTipDir{iMov}{iFrame} .* repmat(objDispVecStartEndUnit(iMov,:),[nTips{iMov}(iFrame) 1]),2));%Minimum angle 
            branchTipAngToDispFull{iMov}{iFrame} = branchTipAng{iMov}{iFrame} - repmat(objDispVecStartEndAngle(iMov,:),[nTips{iMov}(iFrame) 1]);%Full four-quadrant angle
            %Wrap the angles to the expected range.
            branchTipAngToDispFull{iMov}{iFrame}(:,1) = wrapToPi(branchTipAngToDispFull{iMov}{iFrame}(:,1));
            branchTipAngToDispFull{iMov}{iFrame}(:,2) = wrapToPi(branchTipAngToDispFull{iMov}{iFrame}(:,2)*2) / 2;%These go from -pi/2 to pi/2 so we scale before and after wrapping                                    
            
            %And get all branch absolute and relative directions
            
            %First we need to orient edges relative to cell center so we
            %use the distance from this point
            [edgeDistances,closestPt,dToClosest] = analyzeSkeletonDistanceFromPoint(currSkel.vertices,currSkel.edges,currSkel.edgePaths,objTracks{iMov}(iFrame,[2 1 3]) ./ MA(iMov).pixelSize_ );
            [branchDir{iMov}{iFrame},branchAng{iMov}{iFrame}] = calcBranchDirections(currSkel.edges,currSkel.edgePaths,size(currSkel.vertices,1),edgeDistances);
            branchAngToDispFull{iMov}{iFrame} = branchAng{iMov}{iFrame} - repmat(objDispVecStartEndAngle(iMov,:),[numel(currSkel.edgePaths) 1]);%Full four-quadrant angle
            branchAngToDispFull{iMov}{iFrame}(:,1) = wrapToPi(branchAngToDispFull{iMov}{iFrame}(:,1));
            branchAngToDispFull{iMov}{iFrame}(:,2) = wrapToPi(branchAngToDispFull{iMov}{iFrame}(:,2)*2) / 2;%These go from -pi/2 to pi/2 so we scale before and after wrapping                                    
        end        
                                
        meanRadPerMovBranchWeighted(iMov) = mean(meanRadPerFrameBranchWeighted{iMov});
        meanRadPerMovPointWeighted(iMov) = mean(meanRadPerFramePointWeighted{iMov});
        meanTipRadPerMovBranchWeighted(iMov) = mean(meanTipRadPerFrameBranchWeighted{iMov});
        meanTipRadPerMovPointWeighted(iMov) = mean(meanTipRadPerFramePointWeighted{iMov});
        
        meanTotalRadPerMovBranchWeighted(iMov) = mean(totalMeanRadPerFrameBranchWeighted{iMov});
        meanTotalRadPerMovPointWeighted(iMov) = mean(totalMeanRadPerFramePointWeighted{iMov});
        meanTotalTipRadPerMovBranchWeighted(iMov) = mean(totalMeanTipRadPerFrameBranchWeighted{iMov});
        meanTotalTipRadPerMovPointWeighted(iMov) = mean(totalMeanTipRadPerFramePointWeighted{iMov});
        
        
        % ---------- Vertex Degree vs. Motion Figures ---------- %
        
        if p.BatchMode
            vdotFigPM(iMov) = figure('Visible','off');
        else
            vdotFigPM(iMov) = figure;
        end
        
        [ax,h1,h2] = plotyy(tData{iMov}(1:end-1),vertDegreeMean{iMov}(1:end-1),tData{iMov}(1:end-1),objVelPerFrame{iMov});
        set(get(ax(1),'Ylabel'),'String','Mean Vertex Degree')
        set(get(ax(2),'Ylabel'),'String','Velocity, nm/s')
        xlabel('Time, Seconds')        
        saveThatShit('mean vertex degree and centermost point displacement over time',currOutDir)
        
        if p.BatchMode
            vdotFigPM(iMov) = figure('Visible','off');
        else
            vdotFigPM(iMov) = figure;
        end
        
        [ax,h1,h2] = plotyy(tData{iMov}(1:end-1),vertBranchPointMean{iMov}(1:end-1),tData{iMov}(1:end-1),objVelPerFrame{iMov});
        set(get(ax(1),'Ylabel'),'String','Mean Branch Point Degree')
        set(get(ax(2),'Ylabel'),'String','Velocity, nm/s')
        xlabel('Time, Seconds')                                
        saveThatShit('mean branch point degree and centermost point displacement over time',currOutDir)
        
        if p.BatchMode
            vdotFigPM(iMov) = figure('Visible','off');
        else
            vdotFigPM(iMov) = figure;
        end
        
        [ax,h1,h2] = plotyy(tData{iMov}(1:end-1),vertBranchPointTotal{iMov}(1:end-1),tData{iMov}(1:end-1),objVelPerFrame{iMov});
        set(get(ax(1),'Ylabel'),'String','Total Branch Point Degree')
        set(get(ax(2),'Ylabel'),'String','Velocity, nm/s')
        xlabel('Time, Seconds')        
        
        saveThatShit('total branch point degree and centermost point displacement over time',currOutDir)
        
        
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
        
        allBranchTipDir{iMov} = vertcat(branchTipDir{iMov}{:});
        avgDir(iMov,:) = nanmean(allBranchTipDir{iMov},1);
        [avgAng(iMov,1),avgAng(iMov,2)] = cart2sph(avgDir(iMov,1),avgDir(iMov,2),avgDir(iMov,3));
        allBranchTipAng{iMov} = vertcat(branchTipAng{iMov}{:});
        allBranchTipAngToDispFull{iMov} = vertcat(branchTipAngToDispFull{iMov}{:});
        allBranchAngToDispFull{iMov} = vertcat(branchAngToDispFull{iMov}{:});
        branchTipAngHist{iMov}(:,2) = histc(vertcat(allBranchTipAng{iMov}(:,2)),phiBins);
        branchTipAngHist{iMov}(:,1) = histc(vertcat(allBranchTipAng{iMov}(:,1)),thetaBins);
        branchTipAngToDispHist{iMov} = histc(vertcat(branchTipAngToDisp{iMov}{:}),angBins);
        branchTipAngToDispFreq{iMov} = branchTipAngToDispHist{iMov} ./ sum(branchTipAngToDispHist{iMov});%And frequency for averaging
        branchTipAngToDispFullHist{iMov}(:,1) = histc(allBranchTipAngToDispFull{iMov}(:,1),thetaBins);
        branchTipAngToDispFullHist{iMov}(:,2) = histc(allBranchTipAngToDispFull{iMov}(:,2),phiBins);
        branchTipAngToDispFullFreq{iMov} = branchTipAngToDispFullHist{iMov} ./ sum(branchTipAngToDispFullHist{iMov}(:,1));
        
        
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
        
        % ------- Branch Angle / Direction Figures ----- %
        
        %tip only relative angle (minimum angle) histogram figure
        currFig = figure(figArgs{:});
        h = rose(vertcat(branchTipAngToDisp{iMov}{:}),angBins);
        title('Branch tip minimum angle relative to net displacement, all frames');
        figName = [currOutDir filesep 'branch tip minimum angle relative to displacement'];
        mfFigureExport(currFig,figName);
        
        %tip only full four-quadrant angle centered by displacement vector
        currFig = fsFigure(.5);
        subplot(1,2,1);
        h = rose(allBranchTipAngToDispFull{iMov}(:,1),numel(thetaBins));
        title('Branch tip theta relative to net displacement, all frames');
        subplot(1,2,2);
        h = rose(allBranchTipAngToDispFull{iMov}(:,2),numel(phiBins));
        title('Branch phi relative to net displacement, all frames');
        figName = [currOutDir filesep 'branch tip full angle relative to displacement'];
        mfFigureExport(currFig,figName);
        
        
        %All branch full four-quadrant angle centered by displacement vector
        currFig = fsFigure(.5);
        subplot(1,2,1);
        h = rose(allBranchAngToDispFull{iMov}(:,1),numel(thetaBins));
        title('Branch tip theta relative to net displacement, all frames');
        subplot(1,2,2);
        h = rose(allBranchAngToDispFull{iMov}(:,2),numel(phiBins));
        title('Branch phi relative to net displacement, all frames');
        figName = [currOutDir filesep 'branch tip full angle relative to displacement'];
        mfFigureExport(currFig,figName);
        
        
        % ------- Velocity and Branch # time series plots ----- %
        
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
%         
%         if nFramesPerMov(iMov) >= maxLag+1
%             tLags = -(maxLag*timeInt(iMov)):timeInt(iMov):(maxLag*timeInt(iMov));
% 
%             [tmp,~,tmpB] = crosscorr(nTips{iMov}(1:end-1)',objVelPerFrame{iMov},maxLag);
%             figure;
%             plot(tLags,tmp)
%             hold on
%             plot(xlim,ones(1,2)*tmpB(1),'--r')
%             plot(xlim,ones(1,2)*tmpB(2),'--r')
%             plot([0 0 ],ylim,'--k')
%             plot(xlim,[0 0 ],'--k')
%             xlabel('Delay, seconds (Positive Meanse Velocity follows branching')
%             ylabel('Cross Correlation')
%             title('Cross-Corr, Radius-Weighted Branch Number and Instantaneous Velocity')
% 
%             figName = [currOutDir filesep 'cross corr velocity and radius weighted branch number'];
%             print(pOpt{:},[figName '.eps']);
%             print(pOptTIFF{:},[figName '.tif']);
%             hgsave([figName '.fig'])    
% 
%             ccPerCellRadWtVel(iMov,:) = tmp;
%             cbPerCellRadWtVel(iMov,:) = tmpB;
%         
%         else
%             ccPerCellRadWtVel(iMov,:) = NaN;
%             cbPerCellRadWtVel(iMov,:) = NaN;
%         end                
%         
        if nFramesPerMov(iMov) >= maxLag+1

            [tmp,~,tmpB] = crosscorr(vertBranchPointTotal{iMov}(1:end-1)',objVelPerFrame{iMov},maxLag);
            figure;
            plot(tLags,tmp(:,1))
            hold on
            plot(xlim,ones(1,2)*tmpB(1),'--r')
            plot(xlim,ones(1,2)*tmpB(2),'--r')
            plot([0 0 ],ylim,'--k')
            plot(xlim,[0 0 ],'--k')
            xlabel('Delay, seconds (Positive means velocity follows branching')
            ylabel('Cross Correlation')
            title('Cross-Corr, Total Branch Point Degree and Instantaneous Velocity')

            figName = [currOutDir filesep 'cross corr velocity and total branch point degree'];
            print(pOpt{:},[figName '.eps']);
            print(pOptTIFF{:},[figName '.tif']);
            hgsave([figName '.fig'])    

            ccPerCellBranchDegTot(iMov,:) = tmp;
            cbPerCellBranchDegTot(iMov,:) = tmpB;
        else            
            ccPerCellBranchDegTot(iMov,:) = NaN;
            cbPerCellBranchDegTot(iMov,:) = NaN;            
        end
        
        if nFramesPerMov(iMov) >= maxLag+1

            [tmp,~,tmpB] = crosscorr(vertDegreeMean{iMov}(1:end-1)',objVelPerFrame{iMov},maxLag);
            figure;
            plot(tLags,tmp(:,1))
            hold on
            plot(xlim,ones(1,2)*tmpB(1),'--r')
            plot(xlim,ones(1,2)*tmpB(2),'--r')
            plot([0 0 ],ylim,'--k')
            plot(xlim,[0 0 ],'--k')
            xlabel('Delay, seconds (Positive means velocity follows branching')
            ylabel('Cross Correlation')
            title('Cross-Corr, Vertex Degree Mean and Instantaneous Velocity')

            figName = [currOutDir filesep 'cross corr velocity and vertex degree mean'];
            print(pOpt{:},[figName '.eps']);
            print(pOptTIFF{:},[figName '.tif']);
            hgsave([figName '.fig'])    

            ccPerCellVertDegMean(iMov,:) = tmp;
            cbPerCellVertDegMean(iMov,:) = tmpB;
        else            
            ccPerCellVertDegMean(iMov,:) = NaN;
            cbPerCellVertDegMean(iMov,:) = NaN;            
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

%% --------- Branch Angle vs. Motion Combined Distributions -------- %%

%--Tip-only full four-quadrant angle centered by displacement vector ---%%
combBranchTipAngToDispFull = vertcat(allBranchTipAngToDispFull{:});
currFig = fsFigure(.5);
subplot(1,2,1);
h = rose(combBranchTipAngToDispFull(:,1),numel(thetaBins));
title({'Branch tip theta relative to net displacement',['n=' num2str(nMovies) ' cells, m= ' num2str(sum(nFramesPerMov)) ' frames']});
subplot(1,2,2);
h = rose(combBranchTipAngToDispFull(:,2),numel(phiBins));
title('Branch tip phi relative to net displacement, all frames, all movies');
figName = [p.OutputDirectory filesep 'branch tip full angle relative to displacement'];
mfFigureExport(currFig,figName);

%--All-branch four-quadrant angle centered by displacement vector ---%%
combBranchAngToDispFull = vertcat(allBranchAngToDispFull{:});
currFig = fsFigure(.5);
subplot(1,2,1);
h = rose(combBranchAngToDispFull(:,1),numel(thetaBins));
title({'Branch theta relative to net displacement',['n=' num2str(nMovies) ' cells, m= ' num2str(sum(nFramesPerMov)) ' frames']});
subplot(1,2,2);
h = rose(combBranchAngToDispFull(:,2),numel(phiBins));
title('Branch phi relative to net displacement, all frames, all movies');
figName = [p.OutputDirectory filesep 'branch full angle relative to displacement'];
mfFigureExport(currFig,figName);

if ~isempty(p.BranchAngleDisplacementThresh)
    
    actuallyMoves = objDispStartEnd >= p.BranchAngleDisplacementThresh;    
    
    if nnz(actuallyMoves) > 0 
        combBranchAngToDispFullMoving = vertcat(allBranchAngToDispFull{actuallyMoves});

        currFig = fsFigure(.5);
        subplot(1,2,1);
        h = rose(combBranchAngToDispFullMoving(:,1),numel(thetaBins));
        title({'Branch theta relative to net displacement',['n=' num2str(nnz(actuallyMoves)) ' cells, m= ' num2str(sum(nFramesPerMov(actuallyMoves))) ' frames']});
        subplot(1,2,2);
        h = rose(combBranchAngToDispFullMoving(:,2),numel(phiBins));
        title({'Branch phi relative to net displacement, all frames, all movies',['(Both panels use cells with minimum displacement of ' num2str(p.BranchAngleDisplacementThresh/1e3) ' microns)']});
        figName = [p.OutputDirectory filesep 'branch full angle relative to displacement moving cells only'];
        mfFigureExport(currFig,figName);
        
    end
end
%--full four-quadrant angle centered by displacement vector ---%%

% meanBranchTipAngToDispFullFreq = mean(cat(3,branchTipAngToDispFullFreq{:}),3);
% 
%  
% currFig = fsFigure(.5);
% subplot(1,2,1);
% h = roseModTmp(meanBranchTipAngToDispFullFreq(:,1),numel(thetaBins));
% title({'Branch tip theta relative to net displacement',['n=' num2str(nMovies) ' cells, m= ' num2str(sum(nFramesPerMov)) ' frames']});
% subplot(1,2,2);
% h = roseModTmp(meanBranchTipAngToDispFullFreq(:,2),numel(phiBins));
% title('Branch tip phi relative to net displacement, all frames, all movies');
% figName = [currOutDir filesep 'branch tip full angle relative to displacement'];
% mfFigureExport(currFig,figName);



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


%Get combined corr, flip the delays to match with khulouds from above
[meanOfPerCellCCradWt,ciOfPerCellCCradWt] = correlationBootstrap(ccPerCellRadWtVel(isCorrTimeInt,end:-1:1)',cbPerCellRadWtVel(isCorrTimeInt,1)');

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
    ['Average of Per-Cell Correlations, n=' num2str(nnz(isCorrTimeInt))]})

figName = [p.OutputDirectory filesep 'instantaneous velocity and weighted mean tip number avg of per-cell corr'];
print(pOpt{:},[figName '.eps']);
print(pOptTIFF{:},[figName '.tif']);
hgsave([figName '.fig'])    

% [meanOfPerCellCCThreshTip,ciOfPerCellCCThreshTip] = correlationBootstrap(ccPerCellThreshTipVel',cbPerCellThreshTipVel(:,1)');
% 
% figure
% plot(-useTime*maxLag:useTime:useTime*maxLag,meanOfPerCellCCThreshTip)
% hold on
% plot(-useTime*maxLag:useTime:useTime*maxLag,ciOfPerCellCCThreshTip(1,:),'--r')
% legend('Correlation','Boostrapped 95% C.I.')
% plot(-useTime*maxLag:useTime:useTime*maxLag,ciOfPerCellCCThreshTip(2,:),'--r')
% plot(xlim,[ 0 0 ],'--k')
% plot([ 0 0 ],ylim,'--k')
% xlabel('Time Delay, s')
% ylabel('Correlation')
% title({'Thresholded Tip Count and Instantaneous Velocity',...
%     ['Average of Per-Cell Correlations, n=' num2str(nMovies)]})
% 
% figName = [p.OutputDirectory filesep 'instantaneous velocity and thresholded tip number avg of per-cell corr'];
% print(pOpt{:},[figName '.eps']);
% print(pOptTIFF{:},[figName '.tif']);
% hgsave([figName '.fig'])    

%Get combined corr, flip the delays to match with khulouds from above
[meanOfPerCellCCbranchDegTot,ciOfPerCellbranchDegTot] = correlationBootstrap(ccPerCellBranchDegTot(isCorrTimeInt,end:-1:1)',cbPerCellBranchDegTot(isCorrTimeInt,1)');

figure
plot(-useTime*maxLag:useTime:useTime*maxLag,meanOfPerCellCCbranchDegTot)
hold on
plot(-useTime*maxLag:useTime:useTime*maxLag,ciOfPerCellbranchDegTot(1,:),'--r')
legend('Correlation','Boostrapped 95% C.I.')
plot(-useTime*maxLag:useTime:useTime*maxLag,ciOfPerCellbranchDegTot(2,:),'--r')
plot(xlim,[ 0 0 ],'--k')
plot([ 0 0 ],ylim,'--k')
xlabel('Time Delay, s')
ylabel('Correlation')
title({'Total Branch Point Degree and Instantaneous Velocity',...
        'Positive delay means branch degree follows velocity',...
    ['Average of Per-Cell Correlations, n=' num2str(nnz(isCorrTimeInt))]})

figName = [p.OutputDirectory filesep 'instantaneous velocity and total branch point degree avg of per-cell corr'];
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

%% ---- Output ----- %%

%too lazy to list all the relevant variables so just dump it all...
save([p.OutputDirectory filesep 'combined analysis.mat'])

