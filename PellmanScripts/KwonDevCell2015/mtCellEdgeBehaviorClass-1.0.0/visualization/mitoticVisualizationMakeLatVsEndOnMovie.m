function [ output_args ] = mitoticVisualizationMakeLatVsEndOnMovie(projList,varargin)
% mitoticVisualizationMakeLatVsEndOnMovie:   


%% Check INPUT
ip = inputParser;

ip.CaseSensitive = false;

ip.addParameter('outputDirectory',[],@(x) ischar(x)); 

ip.addParameter('minDispVect',3,@(x) isscalar(x));  %in pixels

ip.addParameter('resize',true,@(x) islogical(x)); 

ip.addParameter('title',[],@(x) ischar(x));

ip.addParameter('invert',false,@(x) islogical(x));

ip.addParameter('linewidth',2,@(x) isscalar(x)); 

ip.addParameter('rawToo',true,@(x) islogical(x)); 

ip.addParameter('runffmpeg',true,@(x) islogical(x)); 

ip.addParameter('movieName','movie',@(x) ischar(x)); 

ip.parse(varargin{:});

%% 

saveDirOrig = ip.Results.outputDirectory;
for iProj = 1:numel(projList)
    upDir = upDirectory(projList{iProj},3);
    % upDir2 = upDirectory(projList{iProj},2);
    % get the list of images
    listOfImages = searchFiles('.tif',[],[upDir filesep 'images']);
    %sortList = sortUnpaddedList(listOfImages);
    listOfMasks = searchFiles('.tif',[],[projList{iProj} filesep  'masks']);
    % sortListMasks = sortUnpaddedList(listOfMasks);
    
    [ collectDir,subRoiFolder,cond,cellNum,poleNum ] = getCellIDInfo(projList{iProj});
    saveDir = [saveDirOrig filesep cond cellNum '_' poleNum];
    if ~isdir(saveDir)
        mkdir(saveDir);
    end
    
    if isempty(ip.Results.title)
        title = [cond cellNum 'Pole_' poleNum ]; 
        title =strrep(title,'_',' ');
    end
    
    % load all relavent coordinate structures
    s = load([projList{iProj} filesep 'meta' filesep 'projData.mat']);
    projData = s.projData;
    %    ccFreq = projData.stats.freq;
    s2 = load([projList{iProj} filesep 'meta' filesep 'CorticalInfo' filesep 'corticalData.mat']);
    corticalData = s2.corticalData;
    classMat = corticalData.(['OrientVsDispVsDispMTVect_' (num2str(ip.Results.minDispVect)) ]);
    
    
    
    xMatIn_ByFrame = projData.xMatIn_ByFrame;
    yMatIn_ByFrame = projData.yMatIn_ByFrame;
    dispIn_ByFrame = projData.dispByFrame;
    idxLatLogical = corticalData.(['params' num2str(ip.Results.minDispVect) '_60_0pt7_0pt3']).idxLatLogical; % POTENTIAL VARIABLE CAN VARY WITH TEST
   
    xMatOut_ByFrame = projData.xMatOut_ByFrame;
    yMatOut_ByFrame= projData.yMatOut_ByFrame;
    
    %% resize if necessary
    if ip.Results.resize == 1
        mask = double(imread([listOfMasks{1,2} filesep listOfMasks{1,1}]));
        roiYXCell = bwboundaries(mask);
        roiYX = vertcat(roiYXCell{:});
        minY=floor(min(roiYX(:,1)));
        maxY=ceil(max(roiYX(:,1)));
        minX=floor(min(roiYX(:,2)));
        maxX=ceil(max(roiYX(:,2)));
        % pad
        if minY - 20 < 0
        else
            minY  = minY -20;
        end
        minX = minX - 20;
        maxY = maxY +20;
        maxX = maxX +20;
        
        nyResize = maxY-minY+1;
        nxResize= maxX-minX+1;
        [ny,nx] = size(mask);
        temp=zeros(ny,nx);
        temp(minY:maxY,minX:maxX)=1;
        idx=find(temp(:));
        idxAll=[];
        for i=1:1
            idxAll=[idxAll; idx+(i-1)*(ny*nx)];
        end
        
    end % if resize
    
    
    % combine frames
    xMatIn = vertcat(xMatIn_ByFrame{:});
    yMatIn= vertcat(yMatIn_ByFrame{:});
    dispIn = vertcat(dispIn_ByFrame{:});
    xMatOut = vertcat(xMatOut_ByFrame{:});
    yMatOut = vertcat(yMatOut_ByFrame{:});
    
    if ip.Results.resize == 1 % transform coordinates
        xMatIn  = xMatIn - minX +1;
        yMatIn =  yMatIn - minY +1;
        
        xMatOut = xMatOut - minX +1;
        yMatOut = yMatOut - minY +1;
    end
    
    
    %         % filterIn coordinates by type % have to do this before filter by
    %         frame
    xMatInDispFilt = xMatIn(dispIn>0.3,:) ;
    yMatInDispFilt = yMatIn(dispIn>0.3,:);
    xMatOutFilt = xMatOut(dispIn>0.3,:);
    yMatOutFilt = yMatOut(dispIn>0.3,:);
    %framesFilt  = frameNum(dispIn>0.3,:);
    %
    %
    xMatLat = xMatInDispFilt(idxLatLogical,:);
    yMatLat = yMatInDispFilt(idxLatLogical,:);
    classMatLat = classMat(idxLatLogical,:); % already filtered
    % framesLat = framesFilt(idxLatLogical,:);
    xMatOutLat = xMatOutFilt(idxLatLogical,:);
    yMatOutLat = yMatOutFilt(idxLatLogical,:);
    %
    %
    xMatPer = xMatInDispFilt(~idxLatLogical,:);
    yMatPer = yMatInDispFilt(~idxLatLogical,:);
    classMatPer = classMat(~idxLatLogical,:);
    %  framesPer = framesFilt(~idxLatLogical,:);
    xMatOutPer = xMatOutFilt(~idxLatLogical,:);
    yMatOutPer = yMatOutFilt(~idxLatLogical,:);
    s = load([projList{iProj} filesep 'meta' filesep 'projData.mat']);
    projData = s.projData;
    %         ccFreq = projData.stats.freq;
    s2 = load([projList{iProj} filesep 'meta' filesep 'CorticalInfo' filesep 'corticalData.mat']);
    corticalData = s2.corticalData;
    classMat = corticalData.(['OrientVsDispVsDispMTVect_' (num2str(ip.Results.minDispVect)) ]);
    
    xMatDiscard = xMatIn(dispIn<0.3,:);
    yMatDiscard = yMatIn(dispIn<0.3,:);
    % framesDiscard = frameNum(dispIn<0.3,:);
    xMatOutDiscard = xMatOut(dispIn<0.3,:);
    yMatOutDiscard = yMatOut(dispIn<0.3,:);
    
    nFrames = length(listOfImages(:,1));
    cMap = colormap(cool(100));
    iColor = 40;
    %% Start Overlays
    for iFrame = 1: nFrames
        
        % load first frame for figure
        %imgName = [sortList{iFrame,1} filesep  sortList{iFrame,2} num2str(sortList{iFrame,3})  sortList{iFrame,4}];
        % maskName  = [sortListMasks{iFrame,1} filesep  sortListMasks{iFrame,2} num2str(sortListMasks{iFrame,3})  sortListMasks{iFrame,4}];
        imgName = [listOfImages{iFrame,2} filesep listOfImages{iFrame,1}];
        
        img = double(imread(imgName));
        
        if ip.Results.invert == 1
            
            img = -img;
        end
        mask = double(imread([listOfMasks{iFrame,2} filesep listOfMasks{iFrame,1}]));
        [nyOrig,nxOrig] = size(img);
        if ip.Results.resize == 1
            
            img= reshape(img(idxAll),maxY-minY+1,maxX-minX+1,[]);
            
            mask =  reshape(mask(idxAll),maxY-minY+1,maxX-minX+1,[]);
        end
        
        if ip.Results.rawToo == 1
            img = [img img];
        end
        
        % plot image and mask
        [ny,nx] = size(img);
        setFigure(nx,ny,'off');
        
        imshow(img,[]) ;
        hold on
        roiYX = bwboundaries(mask);
        cellfun(@(x) plot(x(:,2),x(:,1),'color','w'),roiYX);
        hold on
        
        %% LATERAL
        % find lateral in the current frame
        
        idxLatC = (~isnan(xMatLat(:,iFrame)) | ~isnan(xMatOutLat(:,iFrame)));
        xMatLatC = xMatLat(idxLatC,:);
        
        yMatLatC = yMatLat(idxLatC,:);
        xMatOutLatC = xMatOutLat(idxLatC,:);
        yMatOutLatC = yMatOutLat(idxLatC,:);
        if ~isempty(xMatLatC)
            %             % plot the microtubule trajectories by type
            for iTrack = 1:size(xMatLatC,1)
                plot(xMatLatC(iTrack,:),yMatLatC(iTrack,:),'r','LineWidth',ip.Results.linewidth);
                plot(xMatOutLatC(iTrack,:),yMatOutLatC(iTrack,:),'y','LineWidth',ip.Results.linewidth);
                hold on
                %                 if makePlots == 1
                %                     track = xMatLatC(iTrack,:);
                %                     idxLast = find(~isnan(track),1,'last');
                %                     text(xMatLatC(iTrack,idxLast),yMatLatC(iTrack,idxLast),num2str(classMatLat(iTrack,1),2),'Color','y','FontName','Arial');
                %                 end
                hold on
            end
            subTrackLength = length(xMatLatC(:,1));
            % plot connections
            firstPtIdx = arrayfun(@(i) find(~isnan(xMatLatC(i,:)),1,'first'),1:subTrackLength);
            %
            %
            %             % don't plot connectors when
            %             % % extract coords where enters the mask
            xFirst = arrayfun(@(i) xMatLatC(i,firstPtIdx(i)),1:subTrackLength);
            yFirst = arrayfun(@(i) yMatLatC(i,firstPtIdx(i)),1:subTrackLength);
            %
            %
            lastPtIdx = arrayfun(@(i) find(~isnan(xMatOutLatC(i,:)),1,'last'),1:subTrackLength);
            xLast = arrayfun(@(i) xMatOutLatC(i,lastPtIdx(i)),1:subTrackLength);
            yLast= arrayfun(@(i) yMatOutLatC(i,lastPtIdx(i)),1:subTrackLength);
            %
            xConnect = [xFirst' xLast'];
            yConnect = [yFirst' yLast'];
            %
            for iTrack = 1:subTrackLength
                plot(xConnect(iTrack,:),yConnect(iTrack,:),'y','LineWidth',ip.Results.linewidth);
            end
            
            arrayfun(@(x) scatter(xMatLatC(x,iFrame),yMatLatC(x,iFrame),50,'r','+'),1:subTrackLength);
            hold on
            if ip.Results.rawToo == 1
                arrayfun(@(x) quiver(xMatLatC(x,iFrame)-8+nxOrig,yMatLatC(x,iFrame),8,0,0,'Color','r'),1:subTrackLength)
            end
            % plot a small arrow where enters.
            
            % get the frames of xMatLatC
            %frames =  arrayfun(@(x)  find(~isnan(xMatLatC),1,'first'),1:subTrackLength);
            % plot an arrow where there is lateral movement
            
            arrayfun(@(x) scatter(xMatOutLatC(x,iFrame),yMatOutLatC(x,iFrame),50,'y','+'),1:subTrackLength);
            
        end % if ~isempty
        % plot connections
        %% End-ON
        idxPerC = (~isnan(xMatPer(:,iFrame)) | ~isnan(xMatOutPer(:,iFrame)));
        xMatPerC = xMatPer(idxPerC,:);
        
        yMatPerC = yMatPer(idxPerC,:);
        xMatOutPerC = xMatOutPer(idxPerC,:);
        yMatOutPerC = yMatOutPer(idxPerC,:);
        if ~isempty(xMatPerC)
            %             % plot the microtubule trajectories by type
            for iTrack = 1:size(xMatPerC,1)
                plot(xMatPerC(iTrack,:),yMatPerC(iTrack,:),'color',cMap(iColor,:),'LineWidth',ip.Results.linewidth);
                plot(xMatOutPerC(iTrack,:),yMatOutPerC(iTrack,:),'y','LineWidth',ip.Results.linewidth);
                hold on
                %                 if makePlots == 1
                %                     track = xMatLatC(iTrack,:);
                %                     idxLast = find(~isnan(track),1,'last');
                %                     text(xMatLatC(iTrack,idxLast),yMatLatC(iTrack,idxLast),num2str(classMatLat(iTrack,1),2),'Color','y','FontName','Arial');
                %                 end
                hold on
            end
            subTrackLength = length(xMatPerC(:,1));
            % plot connections
            firstPtIdx = arrayfun(@(i) find(~isnan(xMatPerC(i,:)),1,'first'),1:subTrackLength);
            %
            %
            %             % don't plot connectors when
            %             % % extract coords where enters the mask
            xFirst = arrayfun(@(i) xMatPerC(i,firstPtIdx(i)),1:subTrackLength);
            yFirst = arrayfun(@(i) yMatPerC(i,firstPtIdx(i)),1:subTrackLength);
            %
            %
            lastPtIdx = arrayfun(@(i) find(~isnan(xMatOutPerC(i,:)),1,'last'),1:subTrackLength);
            xLast = arrayfun(@(i) xMatOutPerC(i,lastPtIdx(i)),1:subTrackLength);
            yLast= arrayfun(@(i) yMatOutPerC(i,lastPtIdx(i)),1:subTrackLength);
            %
            xConnect = [xFirst' xLast'];
            yConnect = [yFirst' yLast'];
            %
            for iTrack = 1:subTrackLength
                plot(xConnect(iTrack,:),yConnect(iTrack,:),'y','LineWidth',ip.Results.linewidth);
            end
            
            arrayfun(@(x) scatter(xMatPerC(x,iFrame),yMatPerC(x,iFrame),50,cMap(iColor,:),'+'),1:subTrackLength);
            
            arrayfun(@(x) scatter(xMatOutPerC(x,iFrame),yMatOutPerC(x,iFrame),50,'y','+'),1:subTrackLength);
            
        end % if ~isempty
        %% plotDiscard == 1
        plotDiscard = 1;
        if plotDiscard == 1
            idxDiscardC = (~isnan(xMatDiscard(:,iFrame)) | ~isnan(xMatOutDiscard(:,iFrame)));
            xMatDiscardC = xMatDiscard(idxDiscardC,:);
            
            yMatDiscardC = yMatDiscard(idxDiscardC,:);
            xMatOutDiscardC = xMatOutDiscard(idxDiscardC,:);
            yMatOutDiscardC = yMatOutDiscard(idxDiscardC,:);
            
            if ~isempty(xMatDiscardC)
                %             % plot the microtubule trajectories by type
                for iTrack = 1:size(xMatDiscardC,1)
                    plot(xMatDiscardC(iTrack,:),yMatDiscardC(iTrack,:),'w','Linewidth',ip.Results.linewidth);
                    plot(xMatOutDiscardC(iTrack,:),yMatOutDiscardC(iTrack,:),'y','Linewidth',ip.Results.linewidth);
                    hold on
                    %                 if makePlots == 1
                    %                     track = xMatLatC(iTrack,:);
                    %                     idxLast = find(~isnan(track),1,'last');
                    %                     text(xMatLatC(iTrack,idxLast),yMatLatC(iTrack,idxLast),num2str(classMatLat(iTrack,1),2),'Color','y','FontName','Arial');
                    %                 end
                    hold on
                end
                subTrackLength = length(xMatOutDiscardC(:,1));
                % plot connections
                firstPtIdx = arrayfun(@(i) find(~isnan(xMatDiscardC(i,:)),1,'first'),1:subTrackLength);
                %
                %
                %             % don't plot connectors when
                %             % % extract coords where enters the mask
                xFirst = arrayfun(@(i) xMatDiscardC(i,firstPtIdx(i)),1:subTrackLength);
                yFirst = arrayfun(@(i) yMatDiscardC(i,firstPtIdx(i)),1:subTrackLength);
                %
                %
                lastPtIdx = arrayfun(@(i) find(~isnan(xMatOutDiscardC(i,:)),1,'last'),1:subTrackLength);
                xLast = arrayfun(@(i) xMatOutDiscardC(i,lastPtIdx(i)),1:subTrackLength);
                yLast= arrayfun(@(i) yMatOutDiscardC(i,lastPtIdx(i)),1:subTrackLength);
                %
                xConnect = [xFirst' xLast'];
                yConnect = [yFirst' yLast'];
                %
                for iTrack = 1:subTrackLength
                    plot(xConnect(iTrack,:),yConnect(iTrack,:),'y','Linewidth',ip.Results.linewidth);
                end
                
                arrayfun(@(x) scatter(xMatDiscardC(x,iFrame),yMatDiscardC(x,iFrame),50,'w','+'),1:subTrackLength);
                
                arrayfun(@(x) scatter(xMatOutDiscardC(x,iFrame),yMatOutDiscardC(x,iFrame),50,'y','+'),1:subTrackLength);    
            end %
        end %if plotDiscard
        
        %% Extra
        
        % plotScaleBar
        text(nx-50 ,ny-20 ,[num2str(iFrame*projData.secPerFrame-projData.secPerFrame) ' s'],'Color','w','FontName','Arial','FontSize',10);
        pixels = round(5/.108);
        plotScaleBar(pixels,pixels/20,'Label','5 um','Color',[1 1 1]);
        [ny,nx] = size(img);
        text(5,10, title, 'Color','w','FontName','Arial','FontSize',8);
        %                text(nx-100,40,num2str(ccFreq,2),'Color','w','FontName','Arial','FontSize',10);
        if ip.Results.rawToo == 1
            text(nxOrig+5,10,'Red Arrows Show Lateral Moving MTs','Color','r','FontName','Arial','FontSize',9);
        end
        saveas(gcf,[saveDir filesep 'Frame' num2str(iFrame,'%03d') '.eps'],'psc2');
        saveas(gcf,[saveDir filesep 'Frame' num2str(iFrame,'%03d') '.png']);
        close gcf
    end % iframe
    %%
    if ip.Results.runffmpeg
       fr = '5';
       
        fmt = '%03d';
        execute = ['ffmpeg -r ' fr ' -i ' saveDir filesep 'Frame' fmt '.png' ...
            ' -b 2000k ' saveDir filesep ip.Results.movieName  '.wmv'];
        system(execute);
        
        execute = ['ffmpeg -r 5 -i ' saveDir filesep 'Frame' fmt '.png' ...
            ' -crf 22 -b 20000k' saveDir filesep ip.Results.movieName '.mp4'];
        system(execute);
        
%         execute = ['ffmpeg -y -r ' fr ' -i ' fpath filesep 'Frame' fmt '.png '  '-crf 22  -pix_fmt yuv420p '  '-b 50000k -bt 20000k '  fpath filesep movieName '.mp4'];
        cd(saveDir)
        
        %execute = 'mencoder mf://*.png -mf w=800:h=600:fps=5:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o movie.wmv';
        %execute = 'mencoder mf://*.png -mf fps=5:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell:vbitrate=2160000 -oac copy -o movie.wmv';
        
        system(execute);
    end
    
end %

end %

