function [ backboneInfo,frames2Fix,modified,TSFigs ] = GCAneuriteOrientConsistencyCheck(backboneInfo,varargin)
%%  GCAneuriteOrientConsistencyCheck: (STEP II of GCA Segmentation)
%
% This function inputs the backboneInfo.mat data structure output from
% getNeuriteOrientation.m and tests for consistency among the backbone
% seed masks. It makes two main assumptions
% (1) that the neurite entrance point orientation is relatively
% constant throughout the whole movie
% (2) that the majority of backbone seed masks were of the correct
% orientation, so that any outliers found need to be corrected-
% the next best neurite seed candidate that complies with the majority
% vote is chosen

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%    backboneInfo :
%          'backboneInfo': A structure with fields :
%        .backboneSeedMask : an rxc logical mask where r is the height (ny)
%               of the input img and c is the width (nx)
%               marking the ridge that is most likely the neurite entrance.
%               These candidates will be tested for temporal outliers in
%               GCAneuriteOrientConsistencyCheck.m, and these outliers will be
%               corrected.
%               These masks serve as the 'seed' in the first step of
%               GCAveilStemReconstruct.m
%
%       .coordsEnterNeurite : an 1x2 double
%           giving the xy coordinates corresponding to where the neurite
%           entrance point was detected. If there are temporal outliers
%           they will be refined here
%
%
%       .candSeeds : an rxc logical mask where r is the height (ny) of the
%           input img and c is the width (nx)
%           marking any alterative seed ridges that may be considered as
%           entrance cadidates.
%           These other ridges will be re-considered if the selected ridge
%           was found to be a temporal outlier
%
%        .linkedRidgesFinal: an rxc logical mask where r is the height (ny)
%           of the input img and c is the width (nx) marking cleaned large
%           scale ridge candidate paths that will be used for the
%           final reconstruction in GCAveilStemReconstruct.m
% 
% img (OPTIONAL) : rxc array of image  
%        only necessary if TSFigs == true 
%
% PARAMS
%       ('TSOverlays' -> if true will make three troubleshoot
%       folders with the larger outputDirectory for the given channel
%
%       ('SizeOfConsistencyRestraint' -> scalar)
%       Default = 5 pixels
%       This parameter allows the user to control the amount of noise is
%       allowed for the consistency metric. The function currently identifies
%       outliers by simply dilating a small disk mask around the backbone
%       seed entrance points. backbone seed entrance coordinates that fall
%       outside of the resulting connected component are considered
%       potential outliers to be corrected. Smaller values for this parameter
%       can be used when the true neurite entrance varies only slightly throughout
%       the movie, while larger values are sometimes required for neurites
%       backbones that are very mobile.
%
%      'CheckOrient' -> flag to for the user to be able to check and
%       potentialy modify which cluster to align all of the neurite
%       orientations
%       Default = false
%%
%%Input check
ip = inputParser;

ip.CaseSensitive = false;
ip.KeepUnmatched = true;

ip.addRequired('backboneInfo',@(x) isstruct(x));

% Specific
ip.addParameter('TSOverlays',true,@(x) islogical(x));
ip.addParameter('SizeOfConsistencyRestraint',5,@(x) isscalar(x));
ip.addParameter('CheckOrient',false,@(x) islogical(x));

ip.parse(backboneInfo,varargin{:});
p = ip.Results;
%% Initiate
[ySize,xSize] = size(backboneInfo(1).backboneSeedMask);
nFrames = length(backboneInfo);
% initiate modified flag
modified = false;
figCount = 1; 
%% START
coordIn = arrayfun(@(x) x.coordsEnterNeurite,backboneInfo,'uniformoutput',0);
coordInPix = cellfun(@(i) sub2ind([ySize,xSize],i(2),i(1)), coordIn,'uniformoutput',0);
% put into a mask
dilateMask = zeros([ySize,xSize]);
dilateMask(vertcat(coordInPix{:}))=1;
originalPoints = dilateMask;

%makePlotsSpatialClusters  = 0;
dilateMask = imdilate(dilateMask,strel('disk',p.SizeOfConsistencyRestraint)); % changed to 10
if p.CheckOrient == true
    
    roiYX = bwboundaries(dilateMask);
    cellfun(@(x) plot(x(:,2),x(:,1),'color','b','Linewidth',3),roiYX);
    
    hold on
    spy(originalPoints,'g',10);
    hold on
    
end

CCDil  = bwconncomp(dilateMask);
frames2Fix = [];
if(CCDil.NumObjects>1)
    allPixels = vertcat(coordInPix{:});
    idxPixelsInClust = cellfun(@(x) intersect(allPixels,x),CCDil.PixelIdxList,'uniformoutput',0);
    
    numClust = numel(idxPixelsInClust);
    numTimesInClust = nan(numClust,1);
    for iCluster = 1:numClust
        current = idxPixelsInClust{iCluster};
        numTimesInClust(iCluster) =  sum(arrayfun(@(i) sum(allPixels==i),current)) ;
        %         if makePlotsSpatialClusters == 1
        %             [yPlot,xPlot] = ind2sub([ySize,xSize],current(1));
        %
        %             text(xPlot,yPlot,num2str(numTimesInClust(iCluster)),'color','k');
        %             saveas(gcf,[saveDir filesep 'BackboneInputPoints.fig']);
        %         end
    end
    
    % count per frame
    % numPixInClust = cellfun(@(x) length(intersect(allPixels,x)),CCDil.PixelIdxList);
    % pixelsClustLarge = CCDil.PixelIdxList{numPixInClust==max(numPixInClust)};
    pixelsClustLarge = vertcat(idxPixelsInClust{numTimesInClust==max(numTimesInClust)});
    pixelsNOClust = vertcat(idxPixelsInClust{numTimesInClust~=max(numTimesInClust)});
    
    % Ask if the user wants to confirm the location of the neurite
    % entrance. The user can choose a separate point to align
    if ip.Results.CheckOrient == true
        test = zeros(ySize,xSize);
        test(pixelsClustLarge)=1;
        spy(test,'r',10);
        text(10,30,{'Red Marks the Neurite' ; 'Entrance Point Chosen for Alignment'}, ...
            'color','r');
        % roiYX  = bwboundaries(test);
        % cellfun(@(x) plot(x(:,2),x(:,1),'Linewidth',5,'color','b'),roiYX);
        [reply]   =  questdlg('Is the Neurite Orientation Chosen Correct?');
        if strcmpi(reply,'no');
            
            modified = true;
            
            h=impoint;
            position = wait(h);
            %  CC = bwconncomp(dilateMask);
            idx = sub2ind([ySize,xSize],round(position(2)),round(position(1)));
            idxNewClust = cellfun(@(x) (~isempty(intersect(idx,x))), CCDil.PixelIdxList);
            newClust = zeros(ySize,xSize);
            newClust(CCDil.PixelIdxList{idxNewClust})=1;
            roiYX = bwboundaries(newClust);
            cellfun(@(x) plot(x(:,2),x(:,1),'color','m','Linewidth',5), roiYX);
            
            pixelsClustLarge = idxPixelsInClust{idxNewClust};
            pixelsNOClust = vertcat(idxPixelsInClust{~idxNewClust});
            
        end
    end % check orient
    
    test = mat2cell(pixelsNOClust,ones(length(pixelsNOClust),1));
    frames2Fix =  cellfun(@(x) find(allPixels==x),test,'uniformoutput',0);
    frames2Fix = vertcat(frames2Fix{:});
    
end

frames2Fix = sort(frames2Fix);

%% Start Fixing Frames
if isempty(frames2Fix)
    display('All Neurite Orientations Consistent: No Changes Were Made')
    
else
    display('Fixing Minority Neurite Orientations to the Majority')
    
    %% fix backbone seed by finding the ridge closest to the majority cluster
    % make a mask of just the majority cluster
    maskMClust = zeros(ySize,xSize);
    maskMClust(pixelsClustLarge) = 1;
    % dist transform of that mask
    distTransFromMClust = bwdist(maskMClust);
    
    %% get the majority body
        % collect backbone from all good frame
        goodFrames = setdiff(1:nFrames,frames2Fix);
        backbones = arrayfun(@(x) x.backboneSeedMask,backboneInfo(goodFrames),'uniformoutput',0);
        sumBB = zeros(ySize,xSize);
        for iBB = 1:length(goodFrames)
            sumBB = sumBB+backbones{iBB}; % likely better way than a loop for this but quick fix
        end
        
        % find pixels in structure that were more in more than 5 frames (to remove
        % outliers)
        sumBB(sumBB<=5) =0;
        sumBB(sumBB>0) = 1;% make a logical mask
        pixMajBB = find(sumBB==1) ;
        
      
    
    %%
    for iFrame = 1:length(frames2Fix)
        display(['Fixing Frame ' num2str(frames2Fix(iFrame))])
        % load the candidate ridgeMask (mask after cleaning and linear connectivity)
        backboneInfo(frames2Fix(iFrame)).alignmentMask = sumBB; 
        
        candRidges =  backboneInfo(frames2Fix(iFrame)).linkedRidgesFinal;
        
        %
        origBBMask = backboneInfo(frames2Fix(iFrame)).backboneSeedMask;
        % set the original backbone to 0
        candRidges(origBBMask==1) = 0;
        
        
        
        % get the connected component structure of the mask
        CCNMS = bwconncomp(candRidges);
        % get their labels
        labelsRidges = labelmatrix(CCNMS);
        
        % apply the ridge mask
        distRidgeMat = candRidges.*distTransFromMClust;
        %    distances = distances(:);
        %    distances = distances(distances~=0);
        %    minDist = min(distances);
        % find all ridges within 50 pixels
        candLabels = labelsRidges(distRidgeMat<20&distRidgeMat~=0);
        candLabels = unique(candLabels);
     
        
        % find max overlap between cand ridge and
        overlap = cellfun(@(x) length(intersect(pixMajBB,x)),CCNMS.PixelIdxList(candLabels));
        labelKeep = candLabels(overlap==max(overlap));
        
 
        
        %% try to make flag  so if the wrong side was connected can fix. (as in DOCK01 new data)
        if isempty(labelKeep) % sometimes it connects to the wrong side
            % test the original candidate to see if there is significant overlap
            % sometimes simply have a long candidate attach to wrong side
            overlapOrig = length(intersect(pixMajBB,find(origBBMask==1))) ;
            
            % flag to fix
            
            
            
            %
            if ~isempty(overlapOrig)
                % get other endpoint
                EPCoordsBodyMask = getEndpoints(find(origBBMask==1),[ySize,xSize]);
                idxEPs = sub2ind([ySize,xSize],EPCoordsBodyMask(:,2),EPCoordsBodyMask(:,1));
                enterCoords = backboneInfo(frames2Fix(iFrame)).coordsEnterNeurite;
                enterIdx  = sub2ind([ySize,xSize],enterCoords(:,2),enterCoords(:,1));
                
                idxEPs(idxEPs ==enterIdx) = [];
                pixBackboneNew = find(origBBMask==1) ;
                useOtherSide =1 ;
            end
            %      end
        else
            pixBackboneNew = CCNMS.PixelIdxList{labelKeep};
            useOtherSide =0;
        end
        
        
        %    candRidgeIdx = find(cellfun(@(x) ~isempty(intersect(x,pixelsClustLarge)),CCNMS.PixelIdxList));
        %    toTest = CCNMS.PixelIdxList(candRidgeIdx) ;
        %    sizeRidge = cellfun(@(x) length(x),toTest);
        %    candRidgeIdx = candRidgeIdx(sizeRidge==max(sizeRidge));
        %
        %
        
        
        backboneSeed = zeros(ySize,xSize);
        
        
        
        
        backboneSeed(pixBackboneNew)=1;
        
        
        
        
        dims = [ySize,xSize];
        boundaryMask(1:dims(1),1) =1;
        boundaryMask(1:dims(1),dims(2))=1;
        boundaryMask(1,1:dims(2))= 1;
        boundaryMask(dims(1),1:dims(2)) =1;
        if useOtherSide == 1
            backboneSeed(enterIdx) = 0; % make sure to take out the old seed point
        else
        end
        % idxEnterNeurite = find(backboneSeed ==1 & boundaryMask ==1);
        % [yEnter,xEnter] = ind2sub([ySize,xSize],idxEnterNeurite);
        
        %if isempty(idxEnterNeurite) % need to interpolate to make sure closed contour
        % interpolate between to nearest point on boundary
        % find endpoint of backboneSeed
        % filter by the closest distance to the original cluster
        
        EPsBackbone = getEndpoints(pixBackboneNew,[ySize,xSize]);
        idxEPs = sub2ind([ySize,xSize],EPsBackbone(:,2),EPsBackbone(:,1));
        if useOtherSide == 1; % flag to not use old endpoint
            %idxEPs = sub2ind([ySize,xSize],EPsBackbone(:,2),EPsBackbone(:,1));
            enterCoords = backboneInfo(frames2Fix(iFrame)).coordsEnterNeurite;
            enterIdx  = sub2ind([ySize,xSize],enterCoords(:,2),enterCoords(:,1));
            % Take Out the old entrance point
            idxEPs(idxEPs ==enterIdx) = [];
            %[EPY,EPX] = ind2sub([ySize,xSize],idxEPs);
            %EPsBackbone = [EPX, EPY];
        end
        
        EPsDistFromMClust = distTransFromMClust(idxEPs);
        idxEPs = idxEPs(EPsDistFromMClust == min(EPsDistFromMClust));
        [EPY,EPX] = ind2sub([ySize,xSize],idxEPs);
        EPsBackbone = [EPX, EPY];
        [yBoundary,xBoundary] = find(boundaryMask==1);
        [idxBoundaryClose,dist] = KDTreeBallQuery([xBoundary,yBoundary] ,EPsBackbone,20);
        % find the distances from the boundary
        if isempty(idxBoundaryClose{1})
            [idxBoundaryClose,dist] = KDTreeBallQuery([xBoundary,yBoundary],EPsBackbone,100);
        end
        distAll =  vertcat(dist{:});
        
        % find the minimum distance from boundary
        toSave = find(distAll ==min(distAll));
        if length(toSave) > 1
            EPDistFromMClust = distTransFromMClust(idxEPs(toSave));
            toSave = toSave(EPDistFromMClust==min(EPDistFromMClust));
            %toSave = toSave(1);% used to just take the first
        end
        E = arrayfun(@(i) [repmat(i, [numel(idxBoundaryClose{i}) 1]) idxBoundaryClose{i}], 1:length(EPsBackbone(:,1)), 'UniformOutput', false);
        E = vertcat(E{:});
        linkCoords = bresenham([EPsBackbone(E(toSave,1),1),EPsBackbone(E(toSave,1),2)],[xBoundary(E(toSave,2)),yBoundary(E(toSave,2))]);
        backboneSeed(linkCoords(:,2), linkCoords(:,1) ) = 1;
        %linked = 1;
        
        xEnter = xBoundary(E(toSave,2));
        yEnter = yBoundary(E(toSave,2));
        
        %% Optional Troubleshoot Overlay 
% if ip.Results.TSOverlays == true 
%   [ny,nx] = size(backboneSeed); 
%   
%   TSFigs(figCount).h = setFigure(nx,ny,'on'); 
%   TSFigs(figCount).name = 'Re-Alignment';  
%  
%   if ~isempty(ip.Results.img)
%       imshow(-img,[]); 
%       hold on 
%   end 
%  
%     hold on 
%     % plot the information to which to align
%     spy(sumBB,'m'); 
%     % plot the candidate ridges 
%     spy(candRidges,'b'); 
%   
%     spy(
%     
%     % plot the new backbone choice 
%     spy(backboneSeed,''); 
%     
%     
%     figCount = figCount +1; 
%   
% end 
%         
%        
%         
        
        
        
        
        
        clear linkedCoords E
        % end    % need to close contour
        % save new neurite entrance and the new backboneSeed
        backboneInfo(frames2Fix(iFrame)).backboneSeedMask= backboneSeed;
        backboneInfo(frames2Fix(iFrame)).coordsEnterNeurite = [xEnter,yEnter];
         backboneInfo(frames2Fix(iFrame)).timeStamp = clock;
%         hashTag = gcaArchiveGetGitHashTag;
%         backboneInfo(frames2Fix(iFrame)).hashTag = hashTag;
        % document in a folder
        close gcf
        
        
        
    end % iFrame
    
end  % isempty
end % function


