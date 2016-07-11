function [groupList] = mitoticMakeCellEdgeSubroisAndExtractMTTracks(ML,varargin)
% mitoticGetCellEdgeMTTracks: Makes subRegion X Distance from the cell edge.
%                             Makes bipolar subRois if requested using pole positions
%                             Extracts MT tracks from subRegions.
%                             Currently searches for multiple whole cell masks .
%
% From u-track output :       loads :
%                              MD.outputFilePaths_{1} of Comet Detection
%                              Process ('movieInfo')
%                              MD.outputFilesPaths_{1} of Comet Post
%                              Tracking Procss ('projData')
%
% From mitoticGetCellEdgeViaThresholding.m output: loads:
%                              /MD.outputDirectory_/masks  (These are perFrame masks of the cell edge)
%
% If subRoiType  = 'Bipolar'
% From mitoticPerformManualFunctions.m output : loads:
%                             /MD.outputDirectory_/poleInfo.mat (contains information
%                             corresponding to a the point location of
%                             mitotic poles - currently defined manually)
%
% INPUT:
% ML: (REQUIRED) :       movieList (output from movieSelectorGUI.m)
%                        corresponding to the whole cell analysis
%                        that you want to partition 
%
%
% edgeDist: (PARAM) :    positive scalar
%                        Default: 1
%                        The distance from the cell edge (in um) : determines the
%                        width of the subRoi.
%
% subRoiType: (PARAM) : character ('full','bipolar')
%                       'full' (Default)
%                       'bipolar' :    partition the full subRoi into two
%                                      regions based on axis formed by
%                                      the poles 
% bipolarOverlay: (PARAM) : logical
%                          Default: true
%                          If bipolar 'subRoiType' selected will make and save a
%                          sanity overlay.
% 
% outputDirectory (PARAM) : character or empty
%                           Default : [] 
% 
% 
% OUTPUT:
%
%% Check INPUT
ip = inputParser;

ip.CaseSensitive = false;

ip.addParameter('edgeDist',1,@(x) isnumeric(x)); % in um

ip.addParameter('subRoiType','full',@(x) ischar(x));

ip.addParameter('subRoiFilename',[]);

ip.addParameter('bipolarOverlay',true);

ip.addParameter('outputDirectory', pwd); 

ip.addParameter('groupListFilename',[]); 

ip.parse(varargin{:});

%% Set Up 
if isa(ML,'MovieData')
    projList{1,1} = [ML.movieDataPath_ filesep ML.movieDataFileName_];
else
    projList = ML.movieDataFile_;
end

% Set up Default filename if applicable.
if isempty(ip.Results.subRoiFilename)
    subRoiFilename = ['mitoticSubRois_' num2str(ip.Results.edgeDist) '_um_' ip.Results.subRoiType] ;
else
    subRoiFilename = ip.Results.subRoiFilename;
end

switch ip.Results.subRoiType
    case 'full'
       nSubRois = 1;
    case 'bipolar';
       nSubRois = 2;
    otherwise
        error('UnRecognized subRoiType');
end


%% Start

for iProj = 1:numel(projList)
    % Set Up 
    load(projList{iProj}); 
    
    outDir = [MD.outputDirectory_ filesep 'TrackingPackage' filesep 'mtCellEdgeBehaviorClass']; 
    
    if ~isdir(outDir)
        mkdir(outDir);
    end 
        
    display(['Making ' ip.Results.subRoiType  ' subRoi masks for ' projList{iProj}]); 
    
    maskDir = [MD.outputDirectory_ filesep 'masks']; 
    
    % load masks
    [masks] = mitoticSearchFiles('.tif',[],maskDir,0);
    mask1 = [masks{1,2} filesep masks{1,1}];
    
    
    if exist([MD.outputDirectory_ filesep 'roiMaskSpindle.tif'],'file')==0
        
        mask = logical(imread(mask1));
        
        maskInternal = zeros(size(mask));
    else
        maskInternal = logical(imread([MD.outputDirectory_ filesep 'roiMaskSpindle.tif']));
    end
    
    if ip.Results.bipolarOverlay && strcmpi(ip.Results.subRoiType,'bipolar')
        overlayDir = [outDir filesep subRoiFilename filesep 'subRoiOverlays'];
        if ~isdir(overlayDir);
            mkdir(overlayDir);
        end
    end
    %% Run through all frames and make masks
    
    for iFrame = 1:MD.nFrames_
        
       
        fileNameIm = [MD.getChannelPaths{1} filesep MD.getImageFileNames{1}{iFrame}]; 
        img =  double(imread(fileNameIm));
       
        fileNameMask = [masks{iFrame,2} filesep masks{iFrame,1} ];
       
        maskExternal = logical(imread(fileNameMask));
        maskExtraCell = ~maskExternal;
        
        mask = maskExternal-maskInternal;
        
        if sum(maskInternal(:) ) ~= 0
            % get long centroid
            props = regionprops(maskInternal,'centroid','orientation');
            centroid = props.Centroid;
        end
        
        if strcmpi(ip.Results.subRoiType,'bipolar')
            if exist([MD.outputDirectory_ filesep 'poleInfo.mat','file'])==0;
                s = load([MD.outputDirectory_ filesep 'poleInfo.mat']);
            else
                error('No poleInfo.mat found : subRoiType: Bipolar Requires this information');
            end
            
            point = s.poleInfo.coords;
            poles = s.poleInfo.poleMasksInd;
            mAxis = (point(1,2)-point(2,2))/(point(1,1)-point(2,1));
            angle = atand(mAxis);
            bAxis = point(1,2) - mAxis*point(1,1);
            
            mAxisPer = -1/mAxis;
            anglePer = atand(mAxisPer);
            bAxisPer = centroid(2)-mAxisPer*centroid(1);
            
            % full perpendicular ROIS for Asymmetry
            [imL, imW ] = size(mask);
            x = [1:imW];
            y = mAxis +bAxis;
            
            [~,yAll]=meshgrid(1:imW,1:imL);
            yLineAxis=repmat(mAxis.*[1:imW]+bAxis,[imL,1]);
            yPer = repmat(mAxisPer.*[1:imW]+bAxisPer,[imL,1]);
            
            r11 = yAll<=yPer;
            r12 = yAll>yPer;
            
            % this will be used to find intersection and test for long axis label
            x = [1:imW] ;
            y = mAxis.*x+bAxis;
            x = x(y<imL);
            y = y(y<imL);
            x = x(y>0);
            y= y(y>0);
            pixelTest = sub2ind(size(img),ceil(y),x);
            lineMask = zeros(size(mask));
            lineMask(pixelTest) = 1;
            
            % sort the poles to correspond to one pole or the other.
            % pole test
            wholeMask = imfill(mask);
            test1 = r11& wholeMask;
            test2 = r12& wholeMask;
            idxtest1 =  find(test1==1);
            idxpole1 = find(poles(:,:,1)==1);
            overlap= intersect(idxtest1,idxpole1);
            % switch the regions
            if isempty(overlap)
                r11 = yAll>yPer;
                r12 = yAll<=yPer;
            end
        end % if subRoiType = bipolar
        %% SET-UP CORTICAL SUBREGION FOLDERS
        
        %%%%% MAKE SUBROIS %%%%%
        % find initial distance transform from cell edge (1 um)
        subRoiEdgeMask = imfill(mask);
        distTrans = bwdist(~subRoiEdgeMask);
        
        forIntersectAll = subRoiEdgeMask;
        forIntersectAll(distTrans>1) = 0;
        
        
        distTrans = distTrans.*MD.pixelSize_/1000; %  convert to microns
        subRoiEdgeMask(distTrans>ip.Results.edgeDist )=0;
        
        
        if strcmpi(ip.Results.subRoiType,'Bipolar')
            roiSet(:,:,1) = getLargestCC(subRoiEdgeMask & r11);
            roiSet(:,:,2) = getLargestCC(subRoiEdgeMask & r12);
            forIntersect(:,:,1) = forIntersectAll & r11;
            forIntersect(:,:,2) = forIntersectAll & r12;
        else
            roiSet(:,:,1) = subRoiEdgeMask;
        end
        
        for iSub = 1:size(roiSet,3)
            
            subRoiEdgeDir = [outDir filesep subRoiFilename filesep 'sub_' num2str(iSub)];
            if iFrame == 1 % set up new subRoiDirectory
                if ~isdir(subRoiEdgeDir)
                    mkdir(subRoiEdgeDir)
                end
                mkdir([subRoiEdgeDir filesep 'meta'])
                mkdir([subRoiEdgeDir filesep 'feat'])
                movieInfo1  = [MD.processes_{1,1}.outFilePaths_{1}];
                movieInfo2 = [subRoiEdgeDir filesep 'feat' filesep 'movieInfo.mat'];
                copyfile(movieInfo1,movieInfo2);
                
                if ~isdir([subRoiEdgeDir filesep 'masks'])
                    mkdir([subRoiEdgeDir filesep 'masks'])
                end
                
                if ~isdir([subRoiEdgeDir filesep 'extraCell']);
                    mkdir([subRoiEdgeDir filesep 'extraCell']);
                end
                
            end % iFrame ==1
            % save the mask
            fmt = '%03d';
            imwrite(roiSet(:,:,iSub),[subRoiEdgeDir filesep 'masks' filesep 'roiMask' num2str(iFrame, fmt) '.tif']);
            imwrite(maskExtraCell,[subRoiEdgeDir filesep 'extraCell' filesep 'roiMaskEC' num2str(iFrame,fmt) '.tif']);
            
           
        end % iSub
        
       
        
        if ip.Results.bipolarOverlay && strcmpi(ip.Results.subRoiType,'Bipolar')
            [ny,nx] =size(img); 
            setFigure(nx,ny,'off'); 
            imshow(img,[])
            hold on
            plot(x,y,'r');
          
            plotScaleBar( 10*(1000/MD.pixelSize_),5,'Label','10 um');
            colors = ['b','g'];
            for iSub = 1:2
                roiYX = bwboundaries(roiSet(:,:,iSub));
                roiY = roiYX{1}(:,1);
                roiX = roiYX{1}(:,2);
                plot(roiX,roiY,colors(iSub));
                weightedRoi=bwdist(~roiSet(:,:,iSub));
                [r,c]=find(weightedRoi==max(weightedRoi(:)));
                text(c(1),r(1),num2str(iSub),'color','r','fontsize',14);
                
                %                        calculate p (distance from pole to cortex)
                pixBound = find(forIntersect(:,:,iSub) ==1);
                edgePt = intersect(pixBound,pixelTest);
                while isempty(edgePt)
                    dil = zeros(size(mask));
                    dil(pixBound) = 1;
                    dil =  imdilate(dil,strel('disk',1));
                    pixBound = find(dil==1);
                    edgePt = intersect(pixBound,pixelTest);
                end
                if length(edgePt) >1
                    edgePt = edgePt(1); % just take the first
                end
                [yIntersect(iSub), xIntersect(iSub)] = ind2sub(size(img),edgePt);
                scatter(xIntersect(iSub),yIntersect(iSub),100,colors(iSub),'filled');
                d = sqrt((point(iSub,1)-xIntersect(iSub)).^2 + (point(iSub,2)-yIntersect(iSub)).^2);
                p_iFrame = d.*MD.pixelSize_/1000;
                %
                
                text(xIntersect(iSub)+1,yIntersect(iSub)+1,['P = ' num2str(p_iFrame,2) ' um'],'color','y');
                
                poleCortexDist(iSub,iFrame) = p_iFrame;
                
                
            end % for iSub
            %                    plot the poles
            for iPole= 1:2
                spy(poles(:,:,iPole),colors(iPole),30)
            end
            spindleL = sqrt((point(1,1)-point(2,1)).^2 + (point(1,2)-point(2,2)).^2);
            spindleL = spindleL*MD.pixelSize_/1000;
            text(centroid(1)+1,centroid(2)+1,['L = ' num2str(spindleL,2),' um'],'color','y');
        
            %saveas(gcf,[overlayDir filesep num2str(iFrame,'%03d') '.eps'],'psc2');
            %saveas(gcf,[overlayDir filesep num2str(iFrame,'%03d') '.fig']);
            saveas(gcf,[overlayDir filesep num2str(iFrame,'%03d') '.png']);
            close gcf
        end % if bipolarOverlay
        
        clear roiSet forIntersectAll r11 r12 forIntersect
       
    end % for iFrame
    
    if strcmpi(ip.Results.subRoiType,'Bipolar') && ip.Results.bipolarOverlay
        save([MD.outputDirectory_ filesep 'poleCortexDist.mat'],'poleCortexDist'); % make sure to save the calculation for pole cortex dist for each frame
    end %
    clear  poleCortexDist
    display(['Finished Making ' ip.Results.subRoiType  ' mitotic subRoi masks for ' projList{iProj}]); 
 %% Extract MT Tracks for Individual SubRois
    
    for iSub = 1:nSubRois
          
        subRoiEdgeDir = [outDir filesep subRoiFilename filesep 'sub_' num2str(iSub)];
        
        display({'Start Extracting MT Tracks for :' ; subRoiEdgeDir});
        
        % Extract the Tracks
        display({'Start Extracting Tracks for : ' ; subRoiEdgeDir});
        plusTipSubRoiExtractTracksMITOTICPACKAGE_multMasks(MD,subRoiEdgeDir,'fraction',0.01,...
            'extractType','multMasks1','turnFiguresOn',0,'turnHistOn',0,'bipolarMask',...
            strcmpi(ip.Results.subRoiType,'bipolar'));
        
        display({'Finished Extracting MT Tracks for :' ; subRoiEdgeDir});
    end
  
    
end   % for iProj
end

