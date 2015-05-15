function [reconstruct,filoInfo] = gcaAttachFilopodiaStructuresMain(img, skelIn,bodyMask,scaleMap,maxRes,maxTh,analInfoC,normalC,smoothedEdgeC,iFrame,framesPath)
% gcaAttachFilopodiaStructures: was
% cleanMaskWithBackEst_withInternalFiloClean until 20141022
% Cleans the steerable filter ridge response based on the neurite body
% estimation and calls a number of functions to take the messy first pass
% thresholding to an output to a filoInfo data structure and corresponding finished
% reconstruction mask
%
% INPUT: (should condense I think into one structure);
% img = original image needed as we are thresholding out any response that
%       is located in the background: in future can maybe do this step
%       earlier
% skelIn : the response from the steerable filter-maxNMS-thresholding
%          that needs to be cleaned
% bodyMask: the estimated outline from the body of the neurite to make the
% seed (rename veilStem Mask)
% scaleMap: the scalemap corresponding to the multiscale steerable filter %
% no longer using scales (all one scale at this point - was previoulsy
% using the scales to try to get some estimate of the filopodia thickness.
%
% preConnInt
%
% OUTPUT:
% filoInfo: an nx1 structure where n = the number of filopodia in the image
%           (Here a "filopodia" is defined as the endpoint to a branch
%           point): However filopodia can be associated into "tracking
%           objects" which are filopodia branche groups. This indexing occurs dowstream of
%           the data structure.
%           ie each filo has a conIdx (or connectivity idx) and an
%           associated .type.. so far with the nested structures I found
%           this format the most amenable to data extraction and manipulation as opposed to
%           clustering this information "upstream" in the dataStructure: though I might end up changing it
%
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PREPARE HIGH CONFIDENCE RIDGE 'SEEDS' FOR SUBSEQUANT ITERATIVE MATCHING STEPS
% Notes: ridge junctions are typically not reliably detected in the NMS and
% if they are (we should re-check the NMS code - it is debatable if they
% should exist at all)- it is often ambigious as to whether these are a
% cross-over, a branch-point, or noise. Therefore we break them here so we
% can appropriately assign these junction pixels to individual filopdodia in the
% subsequent matching steps.

% Initiate the cleaned array.
skelInClean = skelIn;

% Break the junctions
nn = padarrayXT(double(skelInClean~=0), [1 1]);
sumKernel = [1 1 1];
nn = conv2(sumKernel, sumKernel', nn, 'valid');
nn1 = (nn-1) .* (skelInClean~=0);
junctionMask = nn1>2;
skelInClean(junctionMask) =0;
%testing
% % % figure
% % % imshow(-img,[]);
% % % hold on
% % % spy(skelInClean,'r');
% Remove candidate ridges less than 3 pixels (note this might eventually
% be a problem is have small spans of filopodia between small branches
% Can potentially save these in two ways do not remove if segment
% surrounded by two junctions or could potential clean using graph matching
% later
% individual ridges: connected components
CCRidges = bwconncomp(skelInClean,8); % FIRST PLACE WHERE I BEGIN TO FILTER out signal
csize = cellfun(@(c) numel(c), CCRidges.PixelIdxList);
nsmall = sum(csize<=3);
CCRidges.NumObjects = CCRidges.NumObjects-nsmall;
CCRidges.PixelIdxList(csize<=3) = [];

% MASK OF CLEANED RIDGES MINUS ALL JUNCTIONS
prunedMask = labelmatrix(CCRidges)>0;
%%%spy(prunedMask,'b');
% MASK OF EXTERNAL FILOPODIA RIDGE CANDIDATES
filoTips = prunedMask.*~bodyMask;

% VEIL/STEM MASK (NO FILL)
neuriteEdge = bwboundaries(bodyMask);
edgeMask = zeros(size(img));
idx  = cellfun(@(x) sub2ind(size(img),x(:,1),x(:,2)),neuriteEdge,'uniformoutput',0);
idx = vertcat(idx{:});
edgeMask(idx) = 1;

% CREATE THE HIGH CONDIDENCE 'SEED' (IE THOSE CANDIDATES DIRECTLY ATTACHED
% TO THE VEIL STEM ESTIMATION -
% NOTE: This seed will be used for iterative graph matching steps to
% attach by internal and external filopodia riges based on geometry.
filoExtAll = (filoTips|edgeMask);


%MARIA CHECK BEFORE RELEASE - MAKE SURE OPTION FOR INTERNAL IS STABLE AND
%DOESN'T CRASH
%% INTERNAL LINKING OPTION
preConnInt = 1; % MAKE FORMAL INPUT - NEED TO MAKE AN OPTION
if preConnInt == 1
    filoExtSeedForInt = double(getLargestCC(filoExtAll));
    filoExtSeedForInt = filoExtSeedForInt.*~edgeMask;
    
    %% start to prepare internal filo/ridge candidates for matching
    internalFilo = prunedMask.*bodyMask; %
    internalFiloSpur = bwmorph(internalFilo,'spur',2); % This just cleans things up a bit have to be careful not to lose too much
    % info though sometimes the signal can be a bit weak on some of these so
    % you already only have 2-3 pixels to work with anyway.
    
    % get the connencted components of the internal if CC less than 2 pixels
    % filter
    CCInt = bwconncomp(internalFiloSpur);
    csize = cellfun(@(x) length(x),CCInt.PixelIdxList);
    CCInt.PixelIdxList(csize<2)= [];
    CCInt.NumObjects = CCInt.NumObjects - sum(csize<2);
    
    % keep the internal cand endpoints in a cell array that corresponds to EACH CC
    %  internalFilo1EPs = cellfun(@(x) getEndpoints(x,size(img)),CCInt.PixelIdxList,'uniformoutput',0);
    %
    %
    % %   figure;
    % %  scatter(internalFilo1EPs(:,1),internalFilo1EPs(:,2),'g');
    % %
    %  % filter out those CCs with no or more than 2 end points
    % weirdInt = cellfun(@(x) size(x,1)~=2,internalFilo1EPs);
    % internalFilo1EPs = internalFilo1EPs(~weirdInt);
    % CCInt.PixelIdxList(weirdInt) = [];
    % CCInt.NumObjects = CCInt.NumObjects- sum(weirdInt);
    
    % sanity check
    % internalPix = vertcat(CCInt.PixelIdxList{:});
    % testMask = zeros(size(maxTh));
    % testMask(internalPix) = 1;
    %
    % test = vertcat(internalFilo1EPs{:});
    % scatter(test(:,2),test(:,1),10,'g','filled');
    
    
    
    %% start to prepare high  confidence ext seed filo/ridge candidates for matching
    
    % keep the endpoints in a cell array that corresponds to EACH CC for Seed
    CCFiloExtSeedForInt = bwconncomp(filoExtSeedForInt);
    
    seedFilo1EPs = cellfun(@(x) getEndpoints(x,size(img)),CCFiloExtSeedForInt.PixelIdxList,'uniformoutput',0);
    
    % filter out those CCs with no or more than 2 end points
    weirdSeed = cellfun(@(x) size(x,1)~=2,seedFilo1EPs);
    seedFilo1EPs = seedFilo1EPs(~weirdSeed);
    CCFiloExtSeedForInt.PixelIdxList(weirdSeed) = [];
    CCFiloExtSeedForInt.NumObjects = CCFiloExtSeedForInt.NumObjects- sum(weirdSeed);
    
    %% Prepare for internal matching II:  get only the EPs that are closest to the neurite edge boundary
    
    %  % change EPs to pixIdx
    % [intEPIdx] = cellfun(@(x) sub2ind(size(img),x(:,2),x(:,1)),internalFilo1EPs,'uniformoutput',0);
    %
    % [seedEPIdx] = cellfun(@(x) sub2ind(size(img),x(:,2),x(:,1)),seedFilo1EPs,'uniformoutput',0);
    %
    % % get distTrans relative to neuriteBodyEst
    %  distTrans = bwdist(edgeMask);
    %  [distTransIntEPFromEdge] = cellfun(@(x) distTrans(x),intEPIdx,'uniformoutput',0);
    %  [distTransExtEPFromEdge] = cellfun(@(x) distTrans(x),seedEPIdx,'uniformoutput',0);
    %
    %  [idxKeepInt] = cellfun(@(x) find(x==min(x)),distTransIntEPFromEdge,'uniformoutput',0);
    %
    %  % sometimes might have more than one for each so just take the first
    %  %idxKeepInt = cellfun(@(x) x(1),idxKeepInt);
    %
    %
    %
    %  % find those with dist trans that are exactly the same an indication that
    %  % they are paralllel to the Neurite body edge so filter- might be a faster
    %  % way to do in the future but good enough for now.
    %
    %   parLIdx = cellfun(@(x) size(x,1)>1,idxKeepInt); % get the indices of those that are parallel
    %   % filter these out
    %   idxKeepInt(parLIdx') = [];
    %   CCInt.PixelIdxList(parLIdx') = [];
    %   CCInt.NumObjects = CCInt.NumObjects - sum(parLIdx);
    %
    %
    %  % sometimes might have the two endpoint pixels have the very same
    %  % distTrans therefore need to fix this (could actually remove this here)
    %  % idxKeepExtSeed = cellfun(@(x) x(1), idxKeepExtSeed);
    %
    %   [idxKeepExtSeed] = cellfun(@(x) find(x==min(x)),distTransExtEPFromEdge,'uniformoutput',0);
    %
    % %   parLIdxE = cellfun(@(x) size(x,1)>1,idxKeepExtSeed);
    % %   % filter these out
    % %   idxKeepExtSeed(parLIdxE') = [];
    % %   CCFiloExtSeedForInt.PixelIdxList(parLIdxE') = [];
    %
    %
    %  idxKeepInt = vertcat(idxKeepInt{:});
    %  idxKeepExtSeed = vertcat(idxKeepExtSeed{:})';
    %
    % % for now just do a for loop  ** work out better later if can **
    % % get the internalFiloCoordsClosest to the edge
    % for i = 1:length(idxKeepInt)
    %     internalFilo1EPsFinal{i} = internalFilo1EPs{i}(idxKeepInt(i),:);
    % end
    %
    %  internalFilo1EPsFinal = vertcat(internalFilo1EPsFinal{:}); % only the coords closest to the cell edge will
    %  % be considered for reconstruct
    %
    % for i = 1:length(idxKeepExtSeed)
    %     seedFilo1EPsFinal{i} = seedFilo1EPs{i}(idxKeepExtSeed(i),:);
    % end
    %
    %
    %  seedFilo1EPsFinal = vertcat(seedFilo1EPsFinal{:});
    %
    %  % Only one end-point per candidate/seed (the closet to the neurite body)
    %  % will now be considered for matching
    
    %% Cut internal ridge candidates with strong orientation changes.
    
    % flip the dimensions of the pixIndices such that the point closest to the
    % neuriteEdge is always first so can search for large orientation
    % gradients that indicate the ridge is likely following along the neurite
    % body
    
    %  for i = 1:length(idxKeepInt)
    %        if idxKeepInt(i) == 2;
    %            CCInt.PixelIdxList{i} = flipdim(CCInt.PixelIdxList{i},1);
    %        end
    %  end
    %
    
    % get orientations of internal ridge candidates per pixel from the
    % steerable filter output
    orient = cellfun(@(x) rad2deg(maxTh(x)+pi/2),CCInt.PixelIdxList,'uniformoutput',0);
    
    % calc gradient of orientation
    diffOrient = cellfun(@(x) abs(diff(x)),orient,'uniformoutput',0);
    
    % where orientation differences greater than 90 degrees cut segment.
    orientChangePtsCell= cellfun(@(x) find(x>20 & x <170,1,'first'),diffOrient,'uniformoutput',0);
    
    % for each orientation Change cut
    toChangeVect = cellfun(@(x)  ~isempty(x),orientChangePtsCell);
    IDCCToChange = find(toChangeVect);
    troubleShootFigs.internal= 1;
    % SANITY CHECK
    beforeCut = zeros(size(maxTh));
    beforeCut(vertcat(CCInt.PixelIdxList{:})) = 1;
    if troubleShootFigs.internal == 1;
        [ny,nx] = size(img);
        h = setFigure(nx,ny);
        imshow(img,[]) ;
        hold on
        
        spy(beforeCut,'b');
        hold on
        
    end
    
    
    for iChange = 1:length(IDCCToChange)
        IDC = IDCCToChange(iChange);
        % get the pixel indices to change
        pixIdx = CCInt.PixelIdxList{IDC};
        % get where to cut
        cutHere = orientChangePtsCell{IDCCToChange(iChange)};
        pixIdx(cutHere+1) = [];
        
        CCInt.PixelIdxList{IDC} = pixIdx; % change the pixels
        [y,x] = ind2sub(size(maxTh),pixIdx);
        scatter(x,y,10,'r','filled');
        text(x(1),y(1),num2str(IDC),'color','r');
    end
    
    % plot
    % toPlot = cellfun(@(x) ind2sub(size(maxTh),x),CCInt.PixelIdxList,'uniformoutput',0);
    % cellfun(@(x)
    
    afterCut = zeros(size(maxTh));
    afterCut(vertcat(CCInt.PixelIdxList{:})) = 1;
    %
    %
    if troubleShootFigs.internal ==1;
        
        spy(afterCut,'r');
        cutDir = [framesPath filesep 'Troubleshoot_Internal' filesep 'Cut_Curve'] ;
        
        if ~isdir(cutDir)
            mkdir(cutDir)
        end
        
        saveas(gcf,[cutDir filesep 'Cut_Curvature' num2str(iFrame) '.tif']);
        
    end
    
    %% Reget ConnectedComponents now that have broken based on orientation
    clear CCInt csize
    
    CCInt = bwconncomp(afterCut);
    
    % get rid of singletons
    csize = cellfun(@(x) size(x,1), CCInt.PixelIdxList);
    CCInt.PixelIdxList(csize<=2) =  [];
    CCInt.NumObjects = CCInt.NumObjects - sum(csize<=2);
    
    %GetInternal Endpoints
    
    internalFilo1EPs = cellfun(@(x) getEndpoints(x,size(img)),CCInt.PixelIdxList,'uniformoutput',0);
    %  % change EPs to pixIdx
    
    % make sure to clean out noise (ie those fragments without any endpoints)
    idxLogicNoEPs = cellfun(@(x) isempty(x) ,internalFilo1EPs);
    % get rid of those without endpoints
    internalFilo1EPs(idxLogicNoEPs) = [];
    afterCut(vertcat(CCInt.PixelIdxList{idxLogicNoEPs})) = 0; % set these equal to zero in the original mask
    CCInt.PixelIdxList(idxLogicNoEPs) = [];
    CCInt.NumObjects = CCInt.NumObjects - sum(idxLogicNoEPs);
    
    
    
    % calc displacement vectors: will use for matching (test if better than
    % orientation estimations from steerable filter.
    vectInt =  cellfun(@(x) [x(1,1)-x(2,1), x(1,2) - x(2,2)], internalFilo1EPs ,'uniformoutput',0);
    dInt  = cellfun(@(x) sqrt((x(1,1)-x(2,1))^2 + (x(1,2)-x(2,2))^2),internalFilo1EPs,'uniformoutput',0);
    
    vectSeed =  cellfun(@(x) [x(1,1)-x(2,1), x(1,2) - x(2,2)], seedFilo1EPs ,'uniformoutput',0);
    dSeed = cellfun(@(x) sqrt((x(1,1)-x(2,1))^2 + (x(1,2)-x(2,2))^2),seedFilo1EPs,'uniformoutput',0);
    
    
    
    
    [intEPIdx] = cellfun(@(x) sub2ind(size(img),x(:,2),x(:,1)),internalFilo1EPs,'uniformoutput',0);
    %
    [seedEPIdx] = cellfun(@(x) sub2ind(size(img),x(:,2),x(:,1)),seedFilo1EPs,'uniformoutput',0);
    
    % Don't need to put all pixels through matching only the ones nearest to
    % the cell body should be considered (other option would be to set up to
    % have the two end points compete.  If have time could test processing
    % time for both and see if save anything
    
    % get distTrans relative to neuriteBodyEst
    distTrans = bwdist(edgeMask);
    [distTransIntEPFromEdge] = cellfun(@(x) distTrans(x),intEPIdx,'uniformoutput',0);
    [distTransExtEPFromEdge] = cellfun(@(x) distTrans(x),seedEPIdx,'uniformoutput',0);
    
    [idxKeepInt] = cellfun(@(x) find(x==min(x)),distTransIntEPFromEdge,'uniformoutput',0);
    
    %  % sometimes might have more than one for each so just take the first
    %idxKeepInt = cellfun(@(x) x(1),idxKeepInt);
    %
    %
    %
    %  % find those with dist trans that are exactly the same an indication that
    %  % they are paralllel to the Neurite body edge so filter- might be a faster
    %  % way to do in the future but good enough for now.
    %
    if ~isempty(idxKeepInt) % if no candidates
        parLIdx = cellfun(@(x) size(x,1)>1,idxKeepInt); % get the indices of those that are parallel
        %   % filter these out
        idxKeepInt(parLIdx') = [];
        CCInt.PixelIdxList(parLIdx') = [];
        vectInt(parLIdx')=[];
        dInt(parLIdx') = [];
        CCInt.NumObjects = CCInt.NumObjects - sum(parLIdx);
        internalFilo1EPs(parLIdx') = [];
        
        
        %  Do the same for the external filo
        
        %  % sometimes might have the two endpoint pixels have the very same
        %  % distTrans therefore need to fix this (could actually remove this here)
        %  % idxKeepExtSeed = cellfun(@(x) x(1), idxKeepExtSeed);
        %
        [idxKeepExtSeed] = cellfun(@(x) find(x==min(x)),distTransExtEPFromEdge,'uniformoutput',0);
        %
        
        parLIdxE = cellfun(@(x) size(x,1)>1,idxKeepExtSeed);
        % %   % filter these out
        idxKeepExtSeed(parLIdxE') = [];
        CCFiloExtSeedForInt.PixelIdxList(parLIdxE') = [];
        vectSeed(parLIdxE') =[];
        dSeed(parLIdxE') = [];
        seedFilo1EPs(parLIdxE') = [];
        
        % convert the coordinates to put into graph matching from pixInd to xy
        % coords (put into a cell)
        
        % do it for internal candidate coords
        idxKeepInt = vertcat(idxKeepInt{:});
        if ~isempty(idxKeepInt)
            idxKeepExtSeed = vertcat(idxKeepExtSeed{:})';
            %
            % % for now just do a for loop  ** work out better later if can **
            % % get the internalFiloCoordsClosest to the edge
            for i = 1:length(idxKeepInt)
                internalFilo1EPsFinal{i} = internalFilo1EPs{i}(idxKeepInt(i),:); % just take the point closest to the neurite edge
            end
            %
            %cellfun(@(x) (x(1,1)-x(1,2))^2
            internalFilo1EPsFinal = vertcat(internalFilo1EPsFinal{:}); % only the coords closest to the cell edge will
            %
            % do it for external candidate coords
            for i = 1:length(idxKeepExtSeed)
                seedFilo1EPsFinal{i} = seedFilo1EPs{i}(idxKeepExtSeed(i),:);
            end
            %
            %
            seedFilo1EPsFinal = vertcat(seedFilo1EPsFinal{:});
            %
            %  % Only one end-point per candidate/seed (the closet to the neurite body)
            %  % will now be considered for matching
            
            % make labelMatCandidate
            labelMatCandInt1 = labelmatrix(CCInt);
            
            internalFiloSpur =double(labelMatCandInt1>0);
            reconstruct.Int.Seed{1} = filoExtSeedForInt; % internal candidates
            reconstruct.Int.Cand{1} = internalFiloSpur; % interal Candidates
            
            
            
            
            
            %% Troubleshoot figure internal
            troubleShootFigs.internal = 1;
            if troubleShootFigs.internal == 1
                inSpurDir = [framesPath filesep 'TroubleShoot_Internal' filesep 'internalFiloBeforeSpur' ];
                if ~isdir(inSpurDir)
                    mkdir(inSpurDir)
                end
                [ny,nx] = size(img);
                h = setFigure(nx,ny);
                
                imshow(img,[]);
                hold on
                
                spy(internalFilo,'g'); % green is the original filo detection after thresholding
                saveas(gcf,[inSpurDir filesep 'interalFiloBeforeSpur' num2str(iFrame) '.tif']);
                
                spy(internalFiloSpur,'b'); % blue is after spur
                saveas(gcf,[inSpurDir filesep 'internalFiloAfterCleanUp' num2str(iFrame) '.tif']);
                close gcf
                linkDir = [framesPath filesep 'TroubleShoot_Internal' filesep 'GraphConnect'];
                if ~isdir(linkDir)
                    mkdir(linkDir)
                end
                h = setFigure(nx,ny);
                imshow(img,[]);
                hold on
                
                scatter(internalFilo1EPsFinal(:,1),internalFilo1EPsFinal(:,2),20,'b','filled'); % the endpoint to connect
                scatter(seedFilo1EPsFinal(:,1),seedFilo1EPsFinal(:,2),20,'r','filled');
                spy(internalFiloSpur,'b');
                
                
            end
            %%
            % run through
            % maskpostconnect1 should have all the CC filos after first connection
            % between in and out
            [maskPostConnect1,linkMask1,status]  = gcaConnectInternalFilopodia(internalFilo1EPsFinal,seedFilo1EPsFinal,maxTh,filoExtSeedForInt,10,labelMatCandInt1,vectInt,vectSeed,dInt,dSeed);
            % need to get the internal Filoconnect
            
            
            
            if troubleShootFigs.internal == 1;
                saveas(gcf,[linkDir filesep 'internalInputEPs' num2str(iFrame) '.tif']);
                hold on
                if status ==1
                    spy(linkMask1,'y',10); % plot the link mask
                end
                saveas(gcf,[linkDir filesep 'internalInputEPsWithLinks' num2str(iFrame) '.tif']);
            end
            reconstruct.Int.links{1} = linkMask1;
            reconstruct.Int.Seed{2} = maskPostConnect1;
            
            %% 2nd iteration internal reconstruct
            
            % 2013_07_14 note here you need to likewise make sure that have a 1 to 1
            % attachment (ie both ends do not attach) also need to filter small pieces
            % also might want to check intern filopodia for highest diff in the maxTh
            % along the edge
            
            %%% Will need to change this
            % make an internal filo only mask and get those connected to the edge those
            % will be your high confidence "seed"
            %internalAll = (internalFilo|edgeMask);
            %internalSeed = getLargestCC(internalAll);
            %internalCand = internalAll;
            
            % get the internal seed
            % internalSeed = (maskPostConnect1.*~filoTips) | edgeMask;
            % % try to filter this seed by changes of orientation
            % test= bwconncomp(internalSeed);
            % orientations = cellfun(@(x) maxTh(x), test.PixelIdxList,'uniformoutput',0);
            % % cellfun(@(x) diff(x,1)
            % % find where the derivative is greater than x and cut
            %
            %
            % % get new interal candidates
            % internalFiloSpur(find(internalSeed==1)) =0;
            % % at the very least filter based on pixel size
            % ccInternCand = bwconncomp(internalFiloSpur);
            % csize = cellfun(@(x) length(x), ccInternCand.PixelIdxList);
            % ccInternCand.PixelIdxList(csize<3) = [];
            % internalCand = zeros(size(img)); % initiate mask
            % internalCand(vertcat(ccInternCand.PixelIdxList{:})) = 1; % make new mask of cleaned internal candidates
            %
            % % get endpoints for seed.
            % nn = padarrayXT(double(internalSeed~=0), [1 1]);
            % sumKernel = [1 1 1];
            % nn = conv2(sumKernel, sumKernel', nn, 'valid');
            % nn1 = (nn-1) .* (internalSeed~=0);
            % [EPIntY,EPIntX] = ind2sub(size(img),find(nn1==1));
            % internalSeedEPs = [EPIntX EPIntY];
            
            % clear nn nn1
            %
            % % get endpoints for candidate- NOTE need to have endpoints competing here
            % % to allow only one node connect.
            % % orientation but also potentially distance in this case need to be
            % % associated with the cost. also would want to filter linear structures to
            % % the cell edge.
            % % anything within x distance and
            %
            % nn = padarrayXT(double(internalCand~=0), [1 1]);
            % sumKernel = [1 1 1];
            % nn = conv2(sumKernel, sumKernel', nn, 'valid');
            % nn1 = (nn-1) .* (internalCand~=0);
            % [EPCandY,EPCandX] = ind2sub(size(img),find(nn1==1));
            % internalCandsEPs = [EPCandX EPCandY];
            %
            % labelMatCandInt = bwlabel(internalFiloSpur);
            % final output here needs to be the CCFiloObjs
            % need to make the second iteration such that each EP competing which can
            % be a bit of a pain....try for now to just cut out second iteration and
            % modify how I do the response fits walk out farther and take back
            % farther...
            %[internalMaskPostConnect, linksInternal,status] = connectInternalFilo(internalCandsEPs,internalSeedEPs,maxTh,maskPostConnect1,18,labelMatCandInt);
            
            % reconstruct.Int.Cand{2} = labelMatCandInt>0;
            % reconstruct.Int.links{2}  = internalMaskPostConnect;
            % reconstruct.Int.links{2} = linksInternal;
            % reconstruct.Int.end = internalMaskPostConnect;
            %%
            % if troubleFigInt == 1;
            %     if status == 1;
            %         spy(linksInternal,'y',10);
            %     end
            %     saveas(gcf,[saveDir filesep 'troubleshootInternal' num2str(iFrame,fmt) '.tif']);
            % end
        else
            maskPostConnect1 = filoExtSeedForInt;
        end % if idxKeepInt % 20140301 SEE IF YOU CAN CLEAN THIS UP
        
    else
        maskPostConnect1 = filoExtSeedForInt;
        
    end  %isempty(idxKeepInt)
    
    
else % do not perform internal filopodia matching use the original
    maskPostConnect1 = double(getLargestCC(filoExtAll)); % changed 20141026
    
end % if perform internal matching
%% Record information for the troubleshooting reconstruction movie making
filoSkelPreConnectExt = (filoTips |edgeMask);

reconstruct.input = filoExtAll; % don't input the internal filo for the reconstruct
%% DOCUMENT THE FILOPODIA INFORMATION FROM THE HIGH CONFIDENCE 'SEED' -
% and any internal (ie veil embedded) actin bundles matched in the
% previous step if that option was selected.

CCFiloObjs = bwconncomp(maskPostConnect1);

% that attempts to categorize filo
% filter out small filo
csizeTest = cellfun(@(x) length(x),CCFiloObjs.PixelIdxList);
CCFiloObjs.PixelIdxList(csizeTest<3) = [];
CCFiloObjs.NumObjects = CCFiloObjs.NumObjects - sum(csizeTest<3);

[ filoInfo ] = gcaRecordFilopodiaSeedInformation( CCFiloObjs,img,maxRes,maxTh,edgeMask,bodyMask,analInfoC,normalC,smoothedEdgeC); %% NOTE fix input here!!

%% Reconstruct the external filopodia network from the initial seed

%%%% INITIATE THE ITERATIVE WHILE LOOP %%%%

numViableCand =1; % flag to continue with reconstruction interations as there are viable candidates to attach
filoSkelPreConnect = double(filoSkelPreConnectExt); % initial skeleton before linking : includes all candidates
filoSkelPreConnect = bwmorph(filoSkelPreConnect,'spur');
links = zeros(size(img)); % initiate link matrix
reconIter = 1; % initiate recording reconstruction iterations
linksPre = zeros(size(img));
%%%% BEGIN ITERATING THE REATTACHMENT PROCESS %%%%
status = 1;
while numViableCand >0  % stop the reconstruction process when no more candidates that meet a certain criteria are found
    
    
    % make a label matrix that corresponds to the filoInfo data structure
    % above (this will be updated each iteration)
    labelMatSeedFilo = zeros(size(img));
    pixIndicesSeedFilo = vertcat(filoInfo(:).Ext_pixIndicesBack); % only used the pixel indices measured not the one
    % projected forward
    [yCoordsSeed, xCoordsSeed]= ind2sub(size(img),pixIndicesSeedFilo);
    
    xyCoordsSeedFilo = [xCoordsSeed, yCoordsSeed];
    
    %[xyCoordsSeedFilo] = vertcat(filoInfo(:).coordsXY); % this includes
    %projections forward.
    
    yx  = vertcat(neuriteEdge{:});
    xyCoordsNeurite = [yx(:,2),yx(:,1)];
    
    xySeed = [xyCoordsSeedFilo;xyCoordsNeurite];
    % label the neuriteEdge with 1 and 2
    % create own label mat to ensure your labels are what you think they
    % are!
    for iFilo = 1:numel(filoInfo)
        if ~isnan(filoInfo(iFilo).Ext_pixIndicesBack)
            labelMatSeedFilo(filoInfo(iFilo).Ext_pixIndicesBack) = iFilo+2;
        end
    end
    
    if ~isfield(analInfoC.bodyEst,'pixIndThickBody');
        % add
        thickBodyMask = logical(analInfoC.masks.thickBodyMask);
        analInfoC.bodyEst.pixIndThickBody = find(thickBodyMask==1);
        neuriteEdgeMask = analInfoC.masks.neuriteEdge;
        thinBodyMask = neuriteEdgeMask.*~thickBodyMask;
        analInfoC.bodyEst.pixIndThinBody = find(thinBodyMask==1);
    end
    
    labelMatSeedFilo(analInfoC.bodyEst.pixIndThickBody) = 1; % label thick parts of body as 1
    labelMatSeedFilo(analInfoC.bodyEst.pixIndThinBody) = 2; % label thin parts of body as 2
    
    % labelMatSeedFilo(sub2ind(size(img),neuriteEdge{1}(:,1),neuriteEdge{1}(:,2))) = 1; % the edge is labeled 1
    % this is how you will know if it is backbone
    filoMask = labelMatSeedFilo>0; % heres the new filopodia mask
    reconstruct.seedMask{reconIter} = filoMask;   % always record..
    
    
    
    
    % Get the coordinates of the filo+NeuriteBody 'skel': those filo not attached to
    % body should be included
    %  endpoints = double((filoSkelPreConnect.* (conv2(sumKernel, sumKernel', padarrayXT(filoSkelPreConnect, [1 1]), 'valid')-1))==1); % mask version
    %  [ye,xe] = find(endpoints~=0); % coords of endpoints
    %
    CC = bwconncomp(filoSkelPreConnect|links|linksPre); % connect and refilter candidates
    %labels = double(labelmatrix(CC));
    
    numPix = cellfun(@(x) numel(x),CC.PixelIdxList);
    
    % find and prun the unattached candidates
    %get rid of small CCs lower than x number of pixels and your largest CC
    % maybe make the top ten percent of response etc need to get reattached
    
    %NEW 20141026
    
    maskSeed = zeros([ny,nx]);
    maskSeed(vertcat(CC.PixelIdxList{numPix==max(numPix)}))= 1;
    seedFilos = maskSeed.*~bodyMask;
    CCSeeds = bwconncomp(seedFilos);
    respValuesSeed = cellfun(@(x) maxRes(x),CCSeeds.PixelIdxList,'uniformoutput',0);
    meanRespValuesSeed= cellfun(@(x) mean(x), respValuesSeed);
    sizeSeed = cellfun(@(x) length(x),CCSeeds.PixelIdxList);
    
    % % %     figure;
    % % %     scatter(sizeSeed,meanRespValuesSeed,10,'k','filled');
    
    hold on
    
    CC.PixelIdxList(numPix==max(numPix))= []; % filter out the new seed
    CC.NumObjects = CC.NumObjects -1; % note should make there be a criteria for intensity a
    if reconIter ==1
        
        % get the response of the pieces
        respValuesCand =   cellfun(@(x) maxRes(x),CC.PixelIdxList,'uniformoutput',0);
        meanRespValuesCand = cellfun(@(x) mean(x),respValuesCand);
        sizeCand = cellfun(@(x) length(x), CC.PixelIdxList);
        % get the response of the high confidence seeds
        % % %       scatter(sizeCand,meanRespValuesCand,10,'r','filled');
        % % %       figure
        cutoff = prctile(meanRespValuesSeed,5); % SEE if it is this filtering step here. for frame 120
        % % %
        % % %
        % % %
        % % %       weakCandMask = zeros(ny,nx);
        % % %       strongCandMask = zeros(ny,nx);
        % % %
        toExclude = (meanRespValuesCand<cutoff & sizeCand<10) | sizeCand<=2;
        % % %
        % % %       weakCandMask(vertcat(CC.PixelIdxList{toExclude} ))=1;
        % % %       strongCandMask(vertcat(CC.PixelIdxList{~toExclude}))=1 ;
        % % %       imshow(-img,[])
        % % %       hold on
        % % %       spy(weakCandMask,'r',10);
        % % %       spy(strongCandMask,'g',10);
        % % %       spy(maskSeed,'b',10);
        % % %          close gcf
        
        % erase these candidates completely so they will not be considered in
        % future iterations
        filoSkelPreConnect(vertcat(CC.PixelIdxList{toExclude'}))= 0;
        CC.PixelIdxList(toExclude') = [];
        CC.NumObjects = CC.NumObjects -sum(toExclude);%
        % erase from filoSkelPreConnect
        
    end
    
    % keep on iterating until no more viable candidates
    % check the number of objects here
    if (CC.NumObjects == 0 || status ==0)
        
        break % flag to break loop i snot more viable candidates/linkages
    end
    
    
    % as well
    % prune junctions this will help infinitely later as you will know what
    % type of pieces you will be connnected in the reconstruction will only
    % have two end points
    
    % make mask of candiates for reattachment
    candidateMask1 = labelmatrix(CC) >0;
    candidateMask1 = double(candidateMask1);
    
    if reconIter ==1
        reconstruct.CandMaskPreCluster = candidateMask1;
    end
    % 20141026 I dont' think I need below any more if I did my job
    % correctly above - only would need if mistakenly made junctions in the
    
    % even after first filtering step it helps to prune junctions
    %     nn = padarrayXT(double(candidateMask1~=0), [1 1]);
    %     sumKernel = [1 1 1];
    %     nn = conv2(sumKernel, sumKernel', nn, 'valid');
    %     nn1 = (nn-1) .* (candidateMask1~=0);
    %     junctionMask = nn1>2;
    %     candidateMask1(junctionMask==1) = 0;
    CCCandidates = bwconncomp(candidateMask1);
    %     csize = cellfun(@(x) numel(x),CCCandidates.PixelIdxList);
    %     CCCandidates.PixelIdxList(csize<4)= [];
    %     CCCandidates.NumObjects = CCCandidates.NumObjects - sum(csize<4);
    
    
    
    labelMatCanFilo = labelmatrix(CCCandidates);
    candidateMask1 = labelMatCanFilo>0;
    % now the endpoints are indexed according to the CC Candidates so as
    % unite edges can record the info
    
    
    EPCandidateSort = cellfun(@(x) getEndpoints(x,size(img)),CCCandidates.PixelIdxList,'uniformoutput',0);
    
    if reconIter ==1 % only do this initial clustering step for the first iteration... hmmm need to make sure
        % don't loose this information though
        [ candidateMask1,linkMask,EPCandidateSort,labelMatCanFilo,madeLinks] = connectLinearStructures(EPCandidateSort,maxTh,candidateMask1,labelMatCanFilo,[0.95,0.95,0.95],10);%% 10 
        % it might be easiest just to recalculate the new end points from the
        % new labelMatrix
        % find all post cluster labels these will be slightly different
        % than the pre-cluster labels (ie some will be removed as linear
        % structures are merged) 
        postClustLabels = unique(labelMatCanFilo(labelMatCanFilo~=0));
        numLabels = length(postClustLabels);
        % find all the pixels for these labels 
        pixIdxCand = arrayfun(@(i) find(labelMatCanFilo==postClustLabels(i)), 1:numLabels, 'uniformoutput',0);
        % get each of these labels endpoints
        EPCandidateSort = cellfun(@(x) getEndpoints(x,size(maxTh)),pixIdxCand,'uniformoutput',0);
        
        
        
        reconstruct.CandMaskPostCluster = candidateMask1;
        reconstruct.clusterlinks = linkMask;
        linksPre = linkMask;
        % add these points to the mask
        
    end % if reconIt ==1
    
    
    
    %get rid segments that might not have canonical endpoints as these
    %are very likely noise.
    filoSkelPreConnectFiltered = (filoMask | candidateMask1 );
    % be careful with your indexing here (Note 20150326 yours label numbers
    % maybe different than the length of EPCandidateSort
    nonCanonical = cellfun(@(x) length(x(:,1))~=2,EPCandidateSort);
    
    nonCanonicalIdx = find(nonCanonical==1);
    idxRemove = arrayfun(@(x) find(labelMatCanFilo == postClustLabels(nonCanonicalIdx(x))),1:length(nonCanonicalIdx),'uniformoutput',0);
    labelMatCanFilo(vertcat(idxRemove{:})) = 0;
    EPCandidateSort =  EPCandidateSort(~nonCanonical) ;
    
    % get rid of those with no EPs
    nonEmpty = cellfun(@(x) ~isempty(x),EPCandidateSort);
    EPCandidateSort = EPCandidateSort(nonEmpty);
    % if no more viable candidates break the while loop
    if isempty(EPCandidateSort)
        break
    end
    
    
    % have a gate that is set to 10 pixels so far (things beyond that distance will
    % not be considered)
    % start with 10...
    [outputMasks,filoInfo,status] = gcaConnectExternalFilopodia(xySeed,EPCandidateSort,10,labelMatCanFilo,labelMatSeedFilo,filoSkelPreConnectFiltered,filoInfo,maxRes,maxTh,img,normalC,smoothedEdgeC);
    if status == 1 ;
        %           % note filoInfo will be updated and this will be used to remake the seed
        reconstruct.output{reconIter} = outputMasks;
        links = (links|outputMasks.links);
        
        
        
        
    end % if status
    
    reconIter = reconIter+1; % always go and save new "seed" from data structure even if reconstruction ended
    display(num2str(reconIter))
    
    
end % while
% results1stRound = getLargestCC(filoSkelBranchingFilo);
end

function [coords] = getEndpoints(pixIdx,size)

endpoints = zeros(size);
endpoints(pixIdx)=1;
sumKernel = [1 1 1];
% find endpoints of the floating candidates to attach (note in the
% future might want to prune so that the closest end to the
% body is the only one to be re-attatched: this will avoid double connections)
endpoints = double((endpoints.* (conv2(sumKernel, sumKernel', padarrayXT(endpoints, [1 1]), 'valid')-1))==1);
[ye,xe] = find(endpoints~=0);
coords(:,1) = xe;
coords(:,2) = ye;
end
