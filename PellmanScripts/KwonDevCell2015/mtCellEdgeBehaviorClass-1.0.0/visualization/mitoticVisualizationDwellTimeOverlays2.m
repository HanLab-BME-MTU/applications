function [ output_args ] = mitoticVisualizationDwellTimeOverlays(ML,varargin)
% mitoticVisualizationDwellTimeOverlays
% 
% 
% a groupList which includes paths to the cortical subRoi you would like to
% plot

%% Check Input 
ip = inputParser;

ip.CaseSensitive = false;
ip.addParameter('outputDirectoryCollect',pwd,@(x) ischar(x));

% Parameters for the classification cut-offs
ip.addParameter('minDispVect',3,@(x) isnumeric(x));
ip.addParameter('orientCutOff',60,@(x) isnumeric(x));
ip.addParameter('dispInRegion',0.7,@(x) isnumeric(x));
ip.addParameter('dispDiscard', 0.3,@(x) isnumeric(x));

ip.addParameter('bipolar',true,@(x) islogical(x)); 

ip.parse(varargin{:});

%% set up 
minDispVectStr = num2str(ip.Results.minDispVect);
orientStr = num2str(ip.Results.orientCutOff);
dispInRegionStr = num2str(ip.Results.dispInRegion);
dispDiscardStr= num2str(ip.Results.dispDiscard);

% change any decimals to 'pt'
dispInRegionStr = strrep(dispInRegionStr,'.','pt');
dispDiscardStr = strrep(dispDiscardStr,'.','pt');
fieldnameC = ['params' minDispVectStr '_' orientStr '_' dispInRegionStr '_' dispDiscardStr];

% convert to groupList
input{1} = ML; 
groupList = mitoticCreateGroupListFromMovieLists(input);

%% MAKE DWELL REFINEMENT PLOT
btwGrpNames = unique(groupList,'stable');



%Find the Groups
for iGroup = 1:length(btwGrpNames(:,1))
    idxProjList = strcmp(btwGrpNames{iGroup},groupList(:,1));
    projListC = groupList(idxProjList,:);
    
    % load mask, dwell data, xCoords/yCoords store.
    
    % find combined dwell population values to set cut off
    % determine min/median/max for each population or sort.
    %
    
    for iProj = 1:length(projListC(:,1))
        % Load the extracted tracks for the subregion
        load([projListC{iProj,2} filesep 'meta' filesep 'projData'])
        
        % load the MD 
        load(ML.movieDataFile_{iProj}); 
        
       %  Load the Cortical Classifications
        if exist([projListC{iProj,2} filesep 'meta' filesep 'CorticalInfo' filesep 'corticalData.mat'])==0;
            display('no Cortical MTs : Skipping');
        else
            % load the cortical file 
            load([projListC{iProj,2} filesep 'meta' filesep 'CorticalInfo' filesep 'corticalData.mat']);
            idxLatLog = corticalData.(fieldnameC).idxLatLogical;
            
            % load the image information
            load(ML.movieDataFile_{iProj}); 
            
            dwellsAll = projData.dwellAllTracks;
            disp = vertcat(projData.dispByFrame{:});
            
            % extract pause,shrinkage,terminal,undefined...
            dataMat = projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix; % load to find which tracks end in term,pause, or shrinkage
            dataMat = dataMat(dataMat(:,9) ~= 0,:); % filter the dataMat to get rid of pause, shrinkage, and undefined gaps.
            
            % get the indices for each class (indices are for dataMat filtered)
            idxEndTerm = dataMat(:,9) == 1; % collect the terminal growth events
            idxEndPause = dataMat(:,9) == 2;
            idxEndShrink = dataMat(:,9) == 3; % collect
            
            
            % filter by displacement
            dwells = dwellsAll(disp> ip.Results.dispDiscard & ~isnan(disp));
            % ok I know this is a bit weird here but filter by these when I do the
            % classifications in the cortical region... region where
            % non-optimized
            idxEndTerm = idxEndTerm(disp>ip.Results.dispDiscard & ~isnan(disp));
            idxEndPause = idxEndPause(disp>ip.Results.dispDiscard & ~isnan(disp));
            idxEndShrink = idxEndShrink(disp>ip.Results.dispDiscard & ~isnan(disp));
            
            % get only end on dwells
            dwellsEndOnP = dwells(~idxLatLog & idxEndPause);
            dwellsEndOnST = dwells(~idxLatLog & (idxEndTerm | idxEndShrink));
            
            data = dwellsEndOnST;
            
            % GET IDs for coords - I know it is weird but it helps me not have to change the script
            idxBeforeAllFilt = find(disp>ip.Results.dispDiscard & ~isnan(disp));
            idxBeforeAllFilt = idxBeforeAllFilt((~idxLatLog & (idxEndTerm | idxEndShrink)));
            
            idxPer = idxBeforeAllFilt; % this will be the indexing for the coords
            
            % sanity check for indexing
            
            if ~(isequal(dwellsAll(idxPer),dwellsEndOnST));
                display('Maria you fucked the indexing you retard - check it');
            end
            %% START OVERLAY
            maskFiles = mitoticSearchFiles('.tif',[],[projListC{iProj,2} filesep 'masks'],0,'all',1); 
            roiMask = logical(imread([maskFiles{1}]));
            roiMaskEnd = logical(imread(maskFiles{end}));
            roiYX = bwboundaries(roiMask);
            
            % load image 
            img = double(imread([MD.getChannelPaths{1} filesep MD.getImageFileNames{1}{iFrame}])); 
            
            subRoiEdgeDir = projListC{iProj,2};
            
            xCoordsDwell = projData.xCoordDwell;
            yCoordsDwell = projData.yCoordDwell;
            xCoordsDwellOut = projData.xCoordDwellOutAllTracks;
            yCoordsDwellOut = projData.yCoordDwellOutAllTracks;
            
            if isfield(projData,'xCoordDwellOutAllTracks');
                pad = 40;
            
                minY=floor(min(roiYX{1}(:,1)))-pad;
                maxY=ceil(max(roiYX{1}(:,1)))+pad;
                minX=floor(min(roiYX{1}(:,2)))-pad;
                maxX=ceil(max(roiYX{1}(:,2)))+pad;
                roiMaskC = roiMask;
                roiMaskCEnd = roiMaskEnd;
             
                % crop the image based on roi but avoid making data type double
                
                [imL,imW,dim]=size(img);
                temp=zeros(imL,imW);
                temp(minY:maxY,minX:maxX)=1;
                idx=find(temp(:));
                idxAll=[];
                for i=1:dim
                    idxAll=[idxAll; idx+(i-1)*(imL*imW)];
                end
                
                img=reshape(img(idxAll),maxY-minY+1,maxX-minX+1,[]);
                roiMaskC = reshape(roiMaskC(idxAll),maxY-minY+1,maxX-minX+1,[]);
                roiMaskCEnd = reshape(roiMaskCEnd(idxAll),maxY-minY+1,maxX-minX+1,[]);
                roiYXCrop = bwboundaries(roiMaskC);
                
                
                [ny,nx]= size(img);
                setFigure(nx,ny,'on')
               
                imshow(img,[]);
                hold on
                forScaleBar = 10/0.108;
                plotScaleBar(forScaleBar,forScaleBar/20,'Location','SouthWest','Color',[1,1,1]);
                %plot(roiYXCrop{1}(:,2),roiYXCrop{1}(:,1),'w');
                cellfun(@(x) plot(x(:,2),x(:,1),'w'),roiYXCrop);
                roiYXCropEnd = bwboundaries(roiMaskCEnd);
                cellfun(@(x) plot(x(:,2),x(:,1),'w','Linestyle','--'),roiYXCropEnd);
                load([subRoiEdgeDir filesep 'dwellMasks.mat']);
                %dwellMasks  = dwellMasks(:,:,~projData.startOrEnd);
                dwellMasks = dwellMasks(:,:,idxPer);
                xCoordsInPer = xCoordsDwell(idxPer,:);
                yCoordsInPer = yCoordsDwell(idxPer,:);
                %%
                    % plot connections
                    numDwell = length(dwellsEndOnST);
                    
                    firstPtIdxPer =  arrayfun(@(i) find(~isnan(xCoordsInPer(i,:)),1,'first'),1:length(idxPer));
                    xCoordsFirstIn = arrayfun(@(i) xCoordsInPer(i,firstPtIdxPer(i)),1:length(idxPer));
                    yCoordsFirstIn= arrayfun(@(i) yCoordsInPer(i,firstPtIdxPer(i)),1:length(idxPer));
                    
                    xCoordsOutPer = xCoordsDwellOut(idxPer,:);
                    yCoordsOutPer = yCoordsDwellOut(idxPer,:);
                    lastPtIdxPer = arrayfun(@(i) find(~isnan(xCoordsOutPer(i,:)),1,'last'),1:length(idxPer),'uniformoutput',0);
                    toKeep= cellfun(@(x) ~isempty(x),lastPtIdxPer);
                    lastPtIdxPer = vertcat(lastPtIdxPer{toKeep});
                    xCoordsOutPer = xCoordsOutPer(toKeep,:);
                    yCoordsOutPer = yCoordsOutPer(toKeep,:);
                    
                    xCoordsLastOut = arrayfun(@(i) xCoordsOutPer(i,lastPtIdxPer(i)),1:length(lastPtIdxPer));
                    yCoordsLastOut= arrayfun(@(i) yCoordsOutPer(i,lastPtIdxPer(i)),1:length(lastPtIdxPer));
                    xCoordsFirstIn = xCoordsFirstIn(toKeep);
                    yCoordsFirstIn = yCoordsFirstIn(toKeep);
                    
                    connectX = [xCoordsFirstIn' xCoordsLastOut'];
                    connectY = [yCoordsFirstIn' yCoordsLastOut'];
                    
                    connectX(connectX(:)<minX | connectX(:)>maxX)=nan;
                    connectY(connectY(:)<minY | connectY(:)>maxY)=nan;
                    % transform coords to rectangle frame of reference
                    connectX=connectX-minX+1;
                    connectY=connectY-minY+1;
                    
                    xCoordsDwellOut(xCoordsOutPer(:)<minX | xCoordsOutPer(:)>maxX)=nan;
                    yCoordsDwellOut(yCoordsOutPer(:)<minY | yCoordsOutPer(:)>maxY)=nan;
                    % transform coords to rectangle frame of reference
                    xCoordsOutPer=xCoordsOutPer-minX+1;
                    yCoordsOutPer=yCoordsOutPer-minY+1;
                    
                    for iTrack =1 :length(connectX(:,1))
                        plot(connectX(iTrack,:),connectY(iTrack,:),'w');
                    end
                    
                    for i = 1:numDwell - sum(~toKeep);
                        scatter(xCoordsOutPer(i,:),yCoordsOutPer(i,:),1,'w','filled')
                        
                        plot(xCoordsOutPer(i,:),yCoordsOutPer(i,:),'color','w');
                        
                        %     dwellMaskC = dwellMasks(:,:,i);
                        %     dwellMaskC = reshape(dwellMaskC(idxAll),maxY-minY+1,maxX-minX+1,[]);
                        %     toPlotDwellBox =  bwboundaries(dwellMaskC);
                        %
                        
                        clear toPlotDwellBox
                    end
                    
                xCoordsInPer(xCoordsInPer(:)<minX | xCoordsInPer(:)>maxX)=nan;
                yCoordsInPer(yCoordsInPer(:)<minY | yCoordsInPer(:)>maxY)=nan;
                % transform coords to rectangle frame of reference
                xCoordsInPer=xCoordsInPer-minX+1;
                yCoordsInPer=yCoordsInPer-minY+1;
                
                xCoordsInPer = xCoordsInPer';
                yCoordsInPer= yCoordsInPer';
                
                %% Now that collected all the data define the cMap
                cMapLength=128; cMap=jet(cMapLength);
                
                %data = dwell(idxPer);
                dataRange = 4;
                data(data>dataRange)=dataRange;
                mapper=linspace(0,dataRange,cMapLength)';
                
                % get closest colormap index for each feature
                D=createDistanceMatrix(data,mapper);
                [sD,idxCMap]=sort(abs(D),2);
                
                
                % go through each color of heat map
                for k=1:cMapLength
                    plot(xCoordsInPer(:,idxCMap(:,1)==k),yCoordsInPer(:,idxCMap(:,1)==k),'color',cMap(k,:));
                    if ~isempty(xCoordsInPer(:,idxCMap(:,1)==k));
                        x = xCoordsInPer(:,idxCMap(:,1)==k); y = yCoordsInPer(:,idxCMap(:,1)==k);
                        x = x(:); y= y(:);
                        x = x(~isnan(x)); y = y(~isnan(y));
                        if length(x) == length(y)
                            scatter(x,y,1,cMap(k,:),'filled');
                        end
                    end
                    idxThatColor = find(idxCMap(:,1) ==k);
                    if ~isempty(idxThatColor)
                        for iMask = 1:length(idxThatColor)
                            dwellMaskC = dwellMasks(:,:,idxThatColor(iMask)); % create current mask
                            dwellMaskC = reshape(dwellMaskC(idxAll),maxY-minY+1,maxX-minX+1,[]);
                            toPlotDwellBox =  bwboundaries(dwellMaskC);
                            cellfun(@(x) plot (x(:,2),x(:,1),'color',cMap(k,:)),toPlotDwellBox);
                            clear toPlotDwellBox
                        end % for iMask
                    end % for idxThatColor
                    clear x y
                end % for K
                         
                [~,~,poleNum] = upDirectory(subRoiEdgeDir,1);
                [~,cellNum] = upDirectory(subRoiEdgeDir,4,1);
                [~,perturbation ]= upDirectory(subRoiEdgeDir,5);
                [~,day] = upDirectory(subRoiEdgeDir,6);
                finalSave = [ip.Results.outputDirectoryCollect filesep perturbation]; % collect by condition
                if ~isdir(finalSave)
                    mkdir(finalSave)
                end
                
                saveas(gcf,[finalSave filesep 'RefineDwell_' day '_' cellNum '_' poleNum '.tif']);
                saveas(gcf,[finalSave filesep 'RefineDwell_' day '_' cellNum '_' poleNum '.fig']);
                saveas(gcf,[finalSave filesep 'RefineDwell_' day '_' cellNum '_' poleNum '.eps'],'psc2');
                clear img idxAll toPlot
            end % if isfield
            
            close gcf
            
        end% if exist
    end % iProj
end % iGroup

