function [ output_args ] = mitoticVisualizationDwellTimeOverlays(ML,varargin)
% mitoticVisualizationDwellTimeOverlays : wrapper function make overlays of the dwell times.
% as in Figure 6 E of Kwon et al. Dev Cell 2015
% 
% 
% INPUT: 
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
% clims (PARAM) :  double 
%                  Default [0,4] (s) 
%                  min, max value of the colormap : dwells are colored by 
%                  dwell lifetime. 
% 
% % Parameters for classification cut-offs
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
%% Check Input 
ip = inputParser;

ip.CaseSensitive = false;
ip.addParameter('outputDirectoryCollect',pwd,@(x) ischar(x));

% Parameters for the classification cut-offs
ip.addParameter('minDispVect',3,@(x) isnumeric(x));
ip.addParameter('orientCutOff',60,@(x) isnumeric(x));
ip.addParameter('dispInRegion',0.7,@(x) isnumeric(x));
ip.addParameter('dispDiscard', 0.3,@(x) isnumeric(x));
ip.addParameter('orientDiscard',150,@(x) isnumeric(x)); 

% converts to ML to a subRoiGroupList: therefore need to know if it is 
% looking for the bipolar subRois
ip.addParameter('bipolar',true,@(x) islogical(x)); 


% Plotting options
ip.addParameter('plotScaleBar',true); 
ip.addParameter('pixelSizeMic',0.065); % I set this equal to your pixel size
% note it currently assumes it is the same for all the movies. 

ip.addParameter('invertImage',false); 
ip.addParameter('maskColor',[1,1,1]); 

ip.addParameter('colorScaleBar',[1,1,1]); 
ip.addParameter('lengthOfScaleBar',10); % in um

ip.addParameter('cropToMaskSize',false,@(x) islogical(x)); 
ip.addParameter('padCrop',40); 

ip.addParameter('visible','on'); % 

ip.parse(varargin{:});

%% set up 
minDispVectStr = num2str(ip.Results.minDispVect);
orientStr = num2str(ip.Results.orientCutOff);
dispInRegionStr = num2str(ip.Results.dispInRegion);
dispDiscardStr= num2str(ip.Results.dispDiscard);
orientDiscardStr  = num2str(ip.Results.orientDiscard); 

% change any decimals to 'pt'
dispInRegionStr = strrep(dispInRegionStr,'.','pt');
dispDiscardStr = strrep(dispDiscardStr,'.','pt');
fieldnameC = ['params' minDispVectStr '_' orientStr '_' dispInRegionStr '_' dispDiscardStr '_' orientDiscardStr];

% convert to groupList
input{1} = ML; 
groupList = mitoticCreateGroupListFromMovieLists(input,'bipolar',ip.Results.bipolar);

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
       % load the MD 
      
    for iProj = 1:length(projListC(:,1))
        % Load the extracted tracks for the subregion
        load([projListC{iProj,2} filesep 'meta' filesep 'projData'])
        imgFiles = mitoticSearchFiles('.tif',[],projListC{iProj,3},0,'all',1); 
        img = double(imread(imgFiles{1}));
        if ip.Results.invertImage; 
            img = -img; 
        end 
        
   
       %  Load the Cortical Classifications
        if exist([projListC{iProj,2} filesep 'meta' filesep 'CorticalInfo' filesep 'corticalData.mat'])==0;
            display('no Cortical MTs : Skipping');
        else
            % load the cortical file 
            load([projListC{iProj,2} filesep 'meta' filesep 'CorticalInfo' filesep 'corticalData.mat']);
            idxLatLog = corticalData.(fieldnameC).idxLatLogical;
            classMat = corticalData.(['OrientVsDispVsDispMTVect_' num2str(ip.Results.minDispVect)]); 
            orient = classMat(:,1); 
            
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
            dwells = dwellsAll(disp> ip.Results.dispDiscard & orient< ip.Results.orientDiscard);
            % ok I know this is a bit weird here but filter by these when I do the
            % classifications in the cortical region... region where
            % non-optimized
            idxEndTerm = idxEndTerm(disp>ip.Results.dispDiscard & orient < ip.Results.orientDiscard);
            idxEndPause = idxEndPause(disp>ip.Results.dispDiscard & orient < ip.Results.orientDiscard);
            idxEndShrink = idxEndShrink(disp>ip.Results.dispDiscard & orient < ip.Results.orientDiscard);
            
            % get only end on dwells
            dwellsEndOnP = dwells(~idxLatLog & idxEndPause);
            dwellsEndOnST = dwells(~idxLatLog & (idxEndTerm | idxEndShrink));
            
            data = dwellsEndOnST;
            
            % GET IDs for coords - I know it is weird but it helps me not have to change the script
            idxBeforeAllFilt = find(disp>ip.Results.dispDiscard & orient <ip.Results.orientDiscard);
            idxBeforeAllFilt = idxBeforeAllFilt((~idxLatLog & (idxEndTerm | idxEndShrink)));
            
            idxPer = idxBeforeAllFilt; % this will be the indexing for the coords
            
            % sanity check for indexing
            
            if ~(isequal(dwellsAll(idxPer),dwellsEndOnST));
                display('indexing does not match - check it');
            end
            %% START OVERLAY
            maskFiles = mitoticSearchFiles('.tif',[],[projListC{iProj,2} filesep 'masks'],0,'all',1); 
            roiMask = logical(imread([maskFiles{1}]));
            roiMaskEnd = logical(imread(maskFiles{end}));
            roiYX = bwboundaries(roiMask);
            
           
            subRoiEdgeDir = projListC{iProj,2};
            
            xCoordsDwell = projData.xCoordDwell;
            yCoordsDwell = projData.yCoordDwell;
            xCoordsDwellOut = projData.xCoordDwellOutAllTracks;
            yCoordsDwellOut = projData.yCoordDwellOutAllTracks;
            
            if isfield(projData,'xCoordDwellOutAllTracks');
                roiMaskC = roiMask;
                roiMaskCEnd = roiMaskEnd;
                if ip.Results.cropToMaskSize
                    minY=floor(min(roiYX{1}(:,1)))-ip.Results.padCrop;
                    maxY=ceil(max(roiYX{1}(:,1)))+ip.Results.padCrop;
                    minX=floor(min(roiYX{1}(:,2)))-ip.Results.padCrop;
                    maxX=ceil(max(roiYX{1}(:,2)))+ip.Results.padCrop;
                   
                  
                    
                    % crop the image based on roi but avoid making data type double
                    
                    [imL,imW,dim]=size(img);
                    temp=zeros(imL,imW);
                    temp(minY:maxY,minX:maxX)=1;
                    
                    idxAll=find(temp(:));
                    %                 idxAll=[];
                    %                 for i=1:dim
                    %                     idxAll=[idxAll; idx+(i-1)*(imL*imW)];
                    %                 end
                    
                    img=reshape(img(idxAll),maxY-minY+1,maxX-minX+1,[]);
                    roiMaskC = reshape(roiMaskC(idxAll),maxY-minY+1,maxX-minX+1,[]);
                    roiMaskCEnd = reshape(roiMaskCEnd(idxAll),maxY-minY+1,maxX-minX+1,[]);
                end
                roiYXCrop = bwboundaries(roiMaskC);
                roiYXCropEnd = bwboundaries(roiMaskCEnd);
               
                
                [ny,nx]= size(img);
                setFigure(nx,ny,ip.Results.visible)
               
                imshow(img,[]);
                hold on
                
             
                cellfun(@(x) plot(x(:,2),x(:,1),'color',ip.Results.maskColor),roiYXCrop);
              
                cellfun(@(x) plot(x(:,2),x(:,1),'color',ip.Results.maskColor,'Linestyle','--'),roiYXCropEnd);
                load([subRoiEdgeDir filesep 'dwellMasks.mat']);
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
                    if ip.Results.cropToMaskSize
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
                    end
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
                    if ip.Results.cropToMaskSize
                        xCoordsInPer(xCoordsInPer(:)<minX | xCoordsInPer(:)>maxX)=nan;
                        yCoordsInPer(yCoordsInPer(:)<minY | yCoordsInPer(:)>maxY)=nan;
                        % transform coords to rectangle frame of reference
                        xCoordsInPer=xCoordsInPer-minX+1;
                        yCoordsInPer=yCoordsInPer-minY+1;
                    end
                xCoordsInPer = xCoordsInPer';
                yCoordsInPer= yCoordsInPer';
                
                
                %% Now that collected all the data define the cMap
                cMapLength=128; cMap=jet(cMapLength);
                
                
                data(data>ip.Results.clims(2))=ip.Results.clims(2);
                %
                data(data<ip.Results.clims(1)) = ip.Results.clims(1);
                mapper=linspace(ip.Results.clims(1),ip.Results.clims(2),cMapLength)';
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
                            if ip.Results.cropToMaskSize
                                dwellMaskC = reshape(dwellMaskC(idxAll),maxY-minY+1,maxX-minX+1,[]);
                            end
                            toPlotDwellBox =  bwboundaries(dwellMaskC);
                            cellfun(@(x) plot (x(:,2),x(:,1),'color',cMap(k,:)),toPlotDwellBox);
                            clear toPlotDwellBox
                        end % for iMask
                    end % for idxThatColor
                    clear x y
                end % for K
                
                [~,cellNum] = mitoticUpDirectory(subRoiEdgeDir,5,1);
                
                if ip.Results.bipolar
                    [~,~,poleNum] = mitoticUpDirectory(subRoiEdgeDir,1);
                else 
                    poleNum = 'WholeCell'; 
                end
                
                if ip.Results.plotScaleBar
                  
                    forScaleBar = ip.Results.lengthOfScaleBar/(ip.Results.pixelSizeMic);
                    plotScaleBar(forScaleBar,forScaleBar/20,'Location','SouthWest','Color',ip.Results.colorScaleBar);
                end
                if ip.Results.cropToMaskSize
                    crop = 'Cropped'; 
                else 
                    crop = []; 
                end 
                finalSave = [ip.Results.outputDirectoryCollect filesep 'DwellTimeOverlays' crop... 
                    filesep fieldnameC filesep strrep(ML.movieListFileName_,'.mat','')]; % collect by condition
                if ~isdir(finalSave)
                    mkdir(finalSave)
                end
                
                saveas(gcf,[finalSave filesep 'RefineDwell_'  cellNum '_' poleNum '.tif']);
                saveas(gcf,[finalSave filesep 'RefineDwell_'  cellNum '_' poleNum '.fig']);
                saveas(gcf,[finalSave filesep 'RefineDwell_'  cellNum '_' poleNum '.eps'],'psc2');
                clear img idxAll toPlot
            end % if isfield
            
            close gcf
            
        end% if exist
    end % iProj
end % iGroup

