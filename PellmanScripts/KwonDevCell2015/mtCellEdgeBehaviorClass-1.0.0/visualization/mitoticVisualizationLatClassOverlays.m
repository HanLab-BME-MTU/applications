function [ output_args ] = mitoticVisualizationLatClassOverlays(ML,varargin)
% mitoticVisualizationLatClassOverlays : wrapper function to make overlays of the lateral
% classification as in Figure 7 D of of Kwon et al. Dev. Cell. 2015
%
% 
%% INPUT: 
% movieList: (REQUIRED) : movieList object 
% 
% 
% PARAMETERS: 
% 
% Input/Output
% outputDirectoryCollect (PARAM) : character
%                                  Default : pwd 
%                                  The path where the collected overlays
%                                  will be stored. 
%
%
% bipolar (PARM) : logical 
%                  Default : true 
%                  As we have to convert between movieList with the information 
%                  regarding the input channels and a groupList
%                  with the subRois mitoticCreateGroupListFromMovieLists
%                  requires this information 
%
% indTrackOverlay (PARAM) : logical 
%                           Default : false 
%                           make a folder of overlays for each track
%                           saved in the individual project directory 
% 
% Parameters for classification cut-offs
% minDispVect: (PARAM) : scalar
%                        Default : 3 pixels
%                        the min MT Track Length of interest 
%                        (currently in pixels)
%                        see Supp Figure 1 A of Kwon et al. Dev Cell 2015
%
% orientCutOff: (PARAM) : scalar
%                         Default : 60 (in degrees) 
%                         Local Direction MT Growth relative to the local
%                         cell edge normal.  
%                         lateral classifications as described 
%                         in Figure S7A of Kwon et al. Dev Cell 2015
%                         
%
% dispInRegion: (PARAM) : scalar
%                         Default : 0.7 ( in um) 
%                         MT Trajectories with Cort Dist < this value will 
%                         automatically be classified as end-on. 
%                         see Figure S7A of Kwon et al. Dev Cell 2015
%
% dispDiscard : (PARAM) : scalar 
%                         Default : 0.3 (im um)
%                         MT trajectories wth MT Cort Dist < this value discarded.                         
%                         see Figure S7A of Kwon et al. Dev Cell 2015
% 
% orientDiscard : (PARAM) : scalar 
%                           Default : 150 (in Degrees)
%                           MT trajectories with with Local Direction MT
%                           Growth relative to the local cell edge normal
%                           > than this value will be discarded (likely
%                           error)
%
% 
%
% 
%% Check Input 
ip = inputParser;

ip.CaseSensitive = false;

% input output 
ip.addParameter('outputDirectoryCollect',[],@(x) ischar(x)); % the directory where the collected overlays will be stored
ip.addParameter('bipolar',true,@(x) islogical(x)); 

ip.addParameter('indTrackOverlay',false); % make a folder of overlays for each track 


% Parameters for the classification cut-offs
ip.addParameter('minDispVect',3,@(x) isnumeric(x));
ip.addParameter('orientCutOff',60,@(x) isnumeric(x));
ip.addParameter('dispInRegion',0.7,@(x) isnumeric(x));
ip.addParameter('dispDiscard', 0.3,@(x) isnumeric(x));
ip.addParameter('orientDiscard',150, @(x) isnumeric(x)); 


ip.addParameter('plotScaleBar',true); 
ip.addParameter('pixelSizeMic',0.065); % I set this equal to your pixel size
% note it currently assumes it is the same for all the movies. 

ip.addParameter('invertImage',false); 
ip.addParameter('maskColor',[1,1,1]); 

ip.addParameter('colorScaleBar',[1,1,1]); 
ip.addParameter('lengthOfScaleBar',10); % in um

ip.addParameter('visible','on'); % 
 


ip.parse(varargin{:});
%% Set Up 


% OutputDirectoryCollect 
if isempty(ip.Results.outputDirectoryCollect); 
    saveDir = pwd; 
else 
    saveDir =ip.Results.outputDirectoryCollect;
end 

if ~isdir(saveDir)
    mkdir(saveDir)
end 

minDispVectStr = num2str(ip.Results.minDispVect);
orientStr = num2str(ip.Results.orientCutOff);
dispInRegionStr = num2str(ip.Results.dispInRegion);
dispDiscardStr= num2str(ip.Results.dispDiscard);
orientDiscardStr = num2str(ip.Results.orientDiscard); 

% change any decimals to 'pt'
dispInRegionStr = strrep(dispInRegionStr,'.','pt');
dispDiscardStr = strrep(dispDiscardStr,'.','pt');
fieldnameC = ['params' minDispVectStr '_' orientStr '_' dispInRegionStr '_' dispDiscardStr '_' orientDiscardStr];

% convert to groupList
input{1} = ML; 
groupList = mitoticCreateGroupListFromMovieLists(input,'bipolar',ip.Results.bipolar);
names = groupList(:,1);
btwGrpNames = unique(names,'stable');

%% Start 
for iGroup = 1:length(btwGrpNames)
    % Extract ProjList for current group
    
    projIndxLogic=strcmp(btwGrpNames(iGroup),names);
    projListC = groupList(projIndxLogic,:);
    
    for iProj = 1:size(projListC,1)
        %
        
        toLoad = projListC{iProj,2};
        
        % get info regarding the filename for saving
        [~,cellNum] = mitoticUpDirectory(toLoad,5,1);
        
        if ip.Results.bipolar
            [~,~,poleNum] = mitoticUpDirectory(toLoad,1);
        else
            poleNum = 'WholeCell';
        end
                
        
        listOfImages = mitoticSearchFiles('.tif',[],projListC{iProj,3},0); 
        imgToPlot = double(imread([listOfImages{1,2} filesep listOfImages{1,1}]));
      
        
        % LOAD DATA
        
        load([projListC{iProj,2} filesep 'meta' filesep 'projData.mat']);
        
        load([projListC{iProj,2} filesep 'meta' filesep 'CorticalInfo' filesep 'corticalData.mat']);
        
        classMat = corticalData.(['OrientVsDispVsDispMTVect_' (num2str(ip.Results.minDispVect)) ]);
        orient = classMat(:,1); 
        xMatIn_ByFrame = projData.xMatIn_ByFrame;
        yMatIn_ByFrame = projData.yMatIn_ByFrame;
        dispIn_ByFrame = projData.dispByFrame;
        idxLatLogical = corticalData.(fieldnameC).idxLatLogical; % 
        xMatOut_ByFrame = projData.xMatOut_ByFrame;
        yMatOut_ByFrame= projData.yMatOut_ByFrame;
        
        % find those tracks with coordinates in frame
        frames2Analyze = find(cellfun(@(x) ~isempty(x),xMatIn_ByFrame))';
        
        % repeat the matrix such that can label each track with a frame
        % number
        frameNum = arrayfun(@(i) repmat(frames2Analyze(i),size(xMatIn_ByFrame{frames2Analyze(i)},1),1),1:length(frames2Analyze),'uniformoutput',0);
        frameNum = vertcat(frameNum{:}) ;
         
        % combine frames
        xMatIn = vertcat(xMatIn_ByFrame{:});
        yMatIn= vertcat(yMatIn_ByFrame{:});
        dispIn = vertcat(dispIn_ByFrame{:});
        xMatOut = vertcat(xMatOut_ByFrame{:});
        yMatOut = vertcat(yMatOut_ByFrame{:});
        
        
        % filterIn coordinates by type
        xMatInDispFilt = xMatIn(dispIn>ip.Results.dispDiscard & orient<ip.Results.orientDiscard,:) ;
        yMatInDispFilt = yMatIn(dispIn>ip.Results.dispDiscard & orient < ip.Results.orientDiscard,:);
        xMatOutFilt = xMatOut(dispIn>ip.Results.dispDiscard & orient < ip.Results.orientDiscard,:);
        yMatOutFilt = yMatOut(dispIn>ip.Results.dispDiscard & orient < ip.Results.orientDiscard,:);
        framesFilt  = frameNum(dispIn>ip.Results.dispDiscard & orient < ip.Results.orientDiscard,:);
        
        
        xMatLat = xMatInDispFilt(idxLatLogical,:);
        yMatLat = yMatInDispFilt(idxLatLogical,:);
        classMatLat = classMat(idxLatLogical,:); % already filtered
        framesLat = framesFilt(idxLatLogical,:);
        xMatOutLat = xMatOutFilt(idxLatLogical,:);
        yMatOutLat = yMatOutFilt(idxLatLogical,:);
        
        
        xMatPer = xMatInDispFilt(~idxLatLogical,:);
        yMatPer = yMatInDispFilt(~idxLatLogical,:);
        classMatPer = classMat(~idxLatLogical,:);
        framesPer = framesFilt(~idxLatLogical,:);
        xMatOutPer = xMatOutFilt(~idxLatLogical,:);
        yMatOutPer = yMatOutFilt(~idxLatLogical,:);
        
        xMatDiscard = xMatIn(dispIn<ip.Results.dispDiscard|isnan(dispIn),:);
        yMatDiscard = yMatIn(dispIn<ip.Results.dispDiscard| isnan(dispIn),:);
        framesDiscard = frameNum(dispIn<ip.Results.dispDiscard| isnan(dispIn),:);
        xMatOutDiscard = xMatOut(dispIn<ip.Results.dispDiscard| isnan(dispIn),:);
        yMatOutDiscard = yMatOut(dispIn<ip.Results.dispDiscard| isnan(dispIn),:);
        
        xMatOD = xMatIn(orient>ip.Results.orientDiscard,:); 
        yMatOD = yMatIn(orient>ip.Results.orientDiscard,:);
        framesOD = frameNum(orient>ip.Results.orientDiscard,:); 
        xMatOutOD = xMatOut(orient > ip.Results.orientDiscard,:); 
        yMatOutOD = yMatOut(orient>ip.Results.orientDiscard,:); 
        
        % LOAD MASKS
        maskDir = [toLoad filesep 'masks'] ;
        listOfMasks = mitoticSearchFiles('tif',[],maskDir,0);
        
        mask1 = double(imread([listOfMasks{1,2} filesep listOfMasks{1,1}]));
        roiYX1 = bwboundaries(mask1);
        
        maskend = double(imread([listOfMasks{end,2} filesep listOfMasks{end,1}]));
        roiYXEnd = bwboundaries(maskend);
        
        %% INDIVIDUAL PLOTS
        if ip.Results.indTrackOverlay
            for iFrame = 1:length(frames2Analyze)
                
                frameC = frames2Analyze(iFrame);
                imgC = double(imread([listOfImages{frameC,2} filesep listOfImages{frameC,1}]));
                [ny,nx] = size(imgC) ;
                setFigure(nx,ny,'off');
                imshow(imgC,[]);
                hold on
                % plot currentMask
                
                maskC = double(imread([listOfMasks{frameC,2} filesep listOfMasks{frameC,1}]));
                roiYXC = bwboundaries(maskC);
                
                cellfun(@(x) plot(x(:,2),x(:,1),'color','g'),roiYXC);
                
                % plot just current frame in red or blue
                % get coords in current frame
                
                % lat
                if sum(framesLat==frameC) ~=0
                    numLat = sum(framesLat==frameC);
                    idxLatFrame = find(framesLat==frameC);
                    for iTrack = 1:numLat
                        plot(xMatLat(idxLatFrame(iTrack),:),yMatLat(idxLatFrame(iTrack),:),'r');
                        plot(xMatOutLat(idxLatFrame(iTrack),:),yMatOutLat(idxLatFrame(iTrack),:),'y');
                        %            scatter(xMatLat(idxLatFrame(iTrack),:),yMatLat(idxLatFrame(iFrame),:),20,'r','filled');
                        
                    end
                end
                
                if sum(framesPer == frameC) ~=0
                    numPer = sum(framesPer==frameC);
                    idxPerFrame = find(framesPer==frameC);
                    for iTrack = 1:numPer
                        plot(xMatPer(idxPerFrame(iTrack),:),yMatPer(idxPerFrame(iTrack),:),'b')
                        %scatter(xMatPer(idxPerFrame(iTrack,:),yMatPer(idxPerFrame(iTrack),:),20,'b','filled'));
                        dispC = dispIn_ByFrame{frameC}(iTrack);
                        if dispC < ip.Results.dispInRegion
                            text(nx/2,ny/2,['short disp: ' num2str(dispC, 2)],'color','y','FontSize',14)
                        end
                        plot(xMatOutPer(idxPerFrame(iTrack),:),yMatOutPer(idxPerFrame(iTrack),:),'y');
                    end
                end
                
                if sum(framesDiscard == frameC) ~=0
                    numDiscard = sum(framesDiscard==frameC);
                    idxDiscard = find(framesDiscard==frameC);
                    for iTrack = 1:numDiscard
                        plot(xMatDiscard(idxDiscard(iTrack),:),yMatDiscard(idxDiscard(iTrack),:),'c');
                        % scatter(xMatDiscard(idxDiscard(iTrack),:),yMatDiscard(idxDiscard(iTrack,:)),20,'c','filled');
                        plot(xMatOutDiscard(idxDiscard(iTrack),:),yMatOutDiscard(idxDiscard(iTrack),:),'y');
                    end
                end
                
                if sum(framesOD == frameC) ~=0
                    numOD = sum(framesOD==frameC);
                    idxOD = find(framesOD==frameC);
                    for iTrack = 1:numOD
                        plot(xMatOD(idxOD(iTrack),:),yMatOD(idxOD(iTrack),:),'w');
                        plot(xMatOutOD(idxOD(iTrack),:),yMatOutOD(idxOD(iTrack),:),'y');
                    end
                end
                
                
                outDirI = [projListC{iProj,2} filesep 'meta' filesep 'CorticalInfo' filesep 'TroubleshootClassifications' filesep 'ClassOverlays'... 
                    filesep fieldnameC]; 
                % make a different output directory than this
                            if ~isdir(outDirI);
                                mkdir(outDirI);
                            end
                
                            saveas(gcf,  [outDirI filesep num2str(frameC,'%03d') '.png' ]);
                            saveas(gcf,[outDirI filesep num2str(frameC,'%03d') '.fig']);
                close gcf
            end % sum(frames2analyze) 
        end 
        
        %% FINISH PLOTTING  NOW THAT HAVE PER AND LAT CLASS FINAL: Total
           
            % set figure
            [ny,nx] = size(imgToPlot) ;
            setFigure(nx,ny,'off');
            imshow(imgToPlot,[]);
            hold on
            
            % plot before and after
            cellfun(@(x) plot(x(:,2),x(:,1),'color','g'),roiYX1);
            cellfun(@(x) plot(x(:,2),x(:,1),'--','color','g'),roiYXEnd);
            
            % plot coordinates out
            for iTrack = 1:length(xMatOut(:,1))
                plot(xMatOut(iTrack,:),yMatOut(iTrack,:),'y');
            end
            
            % plot coordinages in different types
            % DISCARDED
            if ~isempty(xMatDiscard)
                for iTrack = 1:size(xMatDiscard,1)
                    plot(xMatDiscard(iTrack,:),yMatDiscard(iTrack,:),'c');
                    hold on
                end
            end % xMatDiscard
            
            % END-ON
            if ~isempty(xMatPer)
                for iTrack = 1:size(xMatPer,1)
                    plot(xMatPer(iTrack,:),yMatPer(iTrack,:),'b');
                    hold on
                end
            end
            
            % LATERAL
            if ~isempty(xMatLat)
                % plot the microtubule trajectories by type
                for iTrack = 1:size(xMatLat,1)
                    plot(xMatLat(iTrack,:),yMatLat(iTrack,:),'r');
                    hold on
                end
            end % if ~isempty
            
            % CONNECTIONS
            subTrackLength = length(xMatIn(:,1));
            %  get idx of the first point of each track in the subregion
            firstPtIdx = arrayfun(@(i) find(~isnan(xMatIn(i,:)),1,'first'),1:subTrackLength);
            
            % don't plot connectors when
            % % extract coords where enters the mask
            xFirst = arrayfun(@(i) xMatIn(i,firstPtIdx(i)),1:subTrackLength);
            yFirst = arrayfun(@(i) yMatIn(i,firstPtIdx(i)),1:subTrackLength);
            
            lastPtIdx = arrayfun(@(i) find(~isnan(xMatOut(i,:)),1,'last'),1:subTrackLength);
            xLast = arrayfun(@(i) xMatOut(i,lastPtIdx(i)),1:subTrackLength);
            yLast= arrayfun(@(i) yMatOut(i,lastPtIdx(i)),1:subTrackLength);
            
            xConnect = [xFirst' xLast'];
            yConnect = [yFirst' yLast'];
            
            for iTrack = 1:subTrackLength
                plot(xConnect(iTrack,:),yConnect(iTrack,:),'y');
            end
           
            % CALC AND PRESENT STATS FOR QUICK REFERENCE
            nLatTot = size(xMatLat,1);
            nPerTot = size(xMatPer,1);
            nDiscard = size(xMatOutDiscard,1);
            percentTrans =     nLatTot/(nLatTot+nPerTot)*100;
            
            text(5,20,['Percent Trans: ' num2str(percentTrans,3)],'color','y');
            text(5,30,['n Side-on: ' num2str(nLatTot) ' Red'],'color','y');
            text(5,40,['n End-on: ' num2str(nPerTot) ' Blue'],'color','y');
            text(5,50,['n Very Low Disp In Region (Discarded):' num2str(nDiscard)],'color','y');
            text(5,60, [ 'n Discard Orient ' num2str(sum(orient>ip.Results.orientDiscard)) ' :not plotted'],'color','y'); 
            text(5,70,['minDispVect: ' num2str(ip.Results.minDispVect) ':Pixels'],'color','y');
            
            condName = strrep(ML.movieListFileName_,'.mat',''); 
            % SAVE  PLOTs
            finalDirectory = [saveDir filesep 'LateralClassOverlays' filesep fieldnameC filesep  condName ];
            
            if ~isdir(finalDirectory)
                mkdir(finalDirectory)
            end
            
            projID = [condName '_Cell_' cellNum  '_' poleNum '_' num2str(ip.Results.minDispVect) '_Pix'];
            saveName = [finalDirectory filesep projID];
            
            
            % save as both a .tif and a .eps
            saveas(gcf,[saveName '.tif']);
            saveas(gcf,[saveName '.eps'],'psc2');
            saveas(gcf,[saveName '.fig']);
            
            close gcf
           
    end % iProj
end % iGroup




