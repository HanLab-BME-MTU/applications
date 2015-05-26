function [ maskPostConnect1,TSFigs,reconstruct ] = gcaReconstructEmbedded(img,maxTh,edgeMask,filoExtSeedForInt,embeddedRidgeCand,varargin)
% gcaReconstructEmbedded: this function connects an external filopodia ridge seed and
% embedded actin bundle ridge candidates based on distance and geometry.  
% Several cleaning steps on the embedded actin bundle ridges are first employed. 
% One of the most important being a cut to orientation. 

%% INPUT: 

% img: (REQUIRED) : RxC double array
%    of image to analyze 
%    where R is the height (ny) and C is the width (nx) of the input image
% 
% maxTh: (REQUIRED) : RxC double array
%    of local orientation of ridge response from steerable filter. 
%    where R is the height (ny) and C is the width (nx) of the input image
%
% edgeMask (REQUIRED) : RxC logical array (binary mask)  
%    of the border pixels of the current veilStem mask 
%    where R is the height (ny) and C is the width
%    (nx) of the original input image
%
% filoExtSeedForInt (REQUIRED) : RxC logical array (binary mask)
%    of the ridge seed associated with  
%    traditional filopodia structures (ie those outside the veil)
%    to which the embedded ridge candidates will be linked
%    where R is the height (ny) and C is the width
%    (nx) of the original input image (Pre-Cleaning Operations)
% 
% embeddedRidgeCand: (REQUIRED) RxC logical array (binary mask) 
%        of the embedded ridge candidates to be linked
%        where R is the height (ny) and C is the width
%        (nx) of the original input image (Pre-Cleaning Operations)
%
%% PARAMS: 
% 'maxRadiusLinkEmbedded' (PARAM) : Scalar 
%          Only embedded ridge candidate end points that are within this max 
%          search radius around each seed ridge endpoint are considered for matching.
%          Default: 10 Pixels
%          See gcaConnectEmbeddedRidgeCandidates.m 
%
%% OUTPUT: 
% maskPostConnect1: RxC logical array (binary mask) 
%        of the final linked mask 
% TSFigs: rx1 structure with fields 
%         where r is the number of troubleshoot figures (to keep adaptive, so 
%         one may add a new TS figure at any point in the code, 
%         the number of troubleshoot figures is not preset. 
%         .h : is the figure handle. 
%         .name : is the character array documenting the name of the
%         figure
%% Parse Input
ip = inputParser;

ip.CaseSensitive = false;
ip.KeepUnmatched = true; 
ip.addRequired('img'); 
ip.addRequired('maxTh'); 
ip.addRequired('filoExtSeedForInt'); 
ip.addRequired('embeddedRidgeCand'); 

ip.addParameter('numPixelsForSpur',2,@(x) isscalar(x)); 

 
ip.addParameter('maxRadiusLinkEmbedded',10,@(x) isscalar(x));
ip.addParameter('TSOverlays',true); 
ip.parse(img,maxTh,filoExtSeedForInt,embeddedRidgeCand,varargin{:});

p = ip.Results;
%% Initiate 
% Initiate figure counter
   countFigs = 1; 
   reconstruct.Int.Seed = []; 
   reconstruct.Int.Cand =[]; 
   reconstruct.Int.links = []; 
%% Clean Embedded Filopodia Candidates

    % This just cleans things up a bit have to be careful not to lose too much
    % info though sometimes the signal can be a bit weak on some of these so
    % you already only have 2-3 pixels to work with anyway.
    embeddedRidgeCandSpur = bwmorph(embeddedRidgeCand,'spur',2);
    
   
    % Remove Singletons
    CCInt = bwconncomp(embeddedRidgeCandSpur);
    csize = cellfun(@(x) length(x),CCInt.PixelIdxList);
    CCInt.PixelIdxList(csize<3)= [];
    CCInt.NumObjects = CCInt.NumObjects - sum(csize<3);% changed to 3 20150526 so to make sure can run through 
    % getEndpoints to get local vectors 
      
%% Clean the Seed 

% take out CCs less than 2 
    CCFiloExtSeedForInt = bwconncomp(filoExtSeedForInt);
    csize = cellfun(@(x) length(x),CCFiloExtSeedForInt.PixelIdxList);
    CCFiloExtSeedForInt.PixelIdxList(csize<3)= [];
    CCFiloExtSeedForInt.NumObjects = CCFiloExtSeedForInt.NumObjects - sum(csize<3);

    
    
    
    [seedFilo1EPs] = cellfun(@(x) getEndpoints(x,size(img),0,1),CCFiloExtSeedForInt.PixelIdxList,'uniformoutput',0);
    
    % filter out those CCs with no or more than 2 end points
    weirdSeed = cellfun(@(x) size(x,1)~=2,seedFilo1EPs);
    seedFilo1EPs = seedFilo1EPs(~weirdSeed);
    
    CCFiloExtSeedForInt.PixelIdxList(weirdSeed) = [];
    CCFiloExtSeedForInt.NumObjects = CCFiloExtSeedForInt.NumObjects- sum(weirdSeed);
    
    % SANITY CHECK
    beforeCut = zeros(size(maxTh));
    beforeCut(vertcat(CCInt.PixelIdxList{:})) = 1;
    
%%  %%  OPTIONAL TS OVERLAY : Before after spur
    if ip.Results.TSOverlays == true;
        [ny,nx] = size(img);
        TSFigs(countFigs).h  = setFigure(nx,ny,'off');
        TSFigs(countFigs).name = 'BeforeAfterCleanStepI';

        imshow(-img,[]) ;
        hold on
        spy(embeddedRidgeCand,'b');
        spy(beforeCut,'g'); 
        text(5,5,'Before Clean Step I','FontSize',10, 'color','b'); 
        text(5,20,'After Clean Step I','FontSize',10,'color','g'); 
        countFigs = countFigs +1; 
    end % 
    
 %%  Cut Embedded Candidates At Points of High Curvature  
 % This step is required because 
 % there are very often embedded ridge candidates that are true candidates 
 % but merge with the ridge response arising at the border of the veil   
 % This helps eliminate this problem.
 
 %% Sort the pixel indices so they are labeled in order from one endpoint 
 % (arbitrary choice) as bwconncomp will not do this completely correctly.
  
 EPsCCsIntForCut = cellfun(@(x) getEndpoints(x,[ny,nx]),CCInt.PixelIdxList,'uniformoutput',0); 
 
 for iCC = 1:numel(CCInt.PixelIdxList) 
     testMask = zeros(ny,nx);  
     testMask(CCInt.PixelIdxList{iCC})=1;
     testMask = logical(testMask);
     startX = EPsCCsIntForCut{iCC}(1,1);
     startY = EPsCCsIntForCut{iCC}(1,2); 
     
     pixIdxBack = nan(length(CCInt.PixelIdxList{iCC}),1);
     transform = bwdistgeodesic(testMask,startX,startY); 
     
      iPix = 0;
        while length(find(transform==iPix)) == 1
            pixIdxBack(iPix+1) = find(transform==iPix); % start at the endpoint
            iPix = iPix +1;
        end
        %pixIdxBack = pixIdxBack(~isnan(pixIdxBack));
     
        CCIntSorted.PixelIdxList{iCC} = pixIdxBack; 
 end 
 CCInt = CCIntSorted; % replace
%%     
  [y,x] = cellfun(@(x) ind2sub([ny,nx],x),CCInt.PixelIdxList,'uniformoutput',0);  
 
   
 
%     if ip.Results.TSOverlays == true 
%         % get the xycoords
         
% 

%% for now make these TS plots as well  
if ip.Results.TSOverlays == true
    %   for iFig = 1:4
    TSFigs(countFigs).h  = setFigure(nx,ny,'on');
    TSFigs(countFigs).name = 'PixelOrder';
    %
    imshow(-img,[]);
    %
    hold on

    for iCC = 1:numel(x)
        x1 = x{iCC};
        y1 = y{iCC};
        %         orient1 = orient{iCC};
        %         diffOrient1 = diffOrient{iCC};
        %
        %         %cmap = jet(length(x1));
        %
        arrayfun(@(i) text(x1(i),y1(i),num2str(i)),1:length(x1));
    end
    countFigs = countFigs +1; 
end
%% Currently simply get orientations of internal ridge candidates per pixel from the
 % steerable filter output. (NOTE MB: investigate better ways to do this
 % before release)
    orient = cellfun(@(x) rad2deg(maxTh(x)+pi/2),CCInt.PixelIdxList,'uniformoutput',0);
     % calc gradient of orientation
    diffOrient = cellfun(@(x) abs(diff(x)),orient,'uniformoutput',0);
%%    
% if ip.Results.TSOverlays == true ; 
%         
% figure; 
% imshow(-img,[]); 
% hold on 
%arrayfun(@(i) plot(x{i},y{i},'color','b'),1:numel(x)); 
         
%         
%         
%         % plotByOrient 
%    cMapLength=128; cMap=jet(cMapLength);
%     orientsAll = vertcat(orient{:});
%                          mapper=linspace(min(orientsAll),max(orientsAll),cMapLength)';
% %                         
% %                         % get closest colormap index for each feature
%                         D=createDistanceMatrix(orientsAll,mapper);
%                         [sD,idxCMap]=sort(abs(D),2);
% %                         
% for iCC = 1:numel(x) % for each piece
%     D = 
%     % plot the orientation 
%     arrayfun(@(i) plot(x
% end 
% end         

%%    
    % originally 20 
    % where orientation differences greater than 90 degrees cut segment.
    orientChangePtsCell= cellfun(@(x) find(x>20 & x <170,1,'first'),diffOrient,'uniformoutput',0);
    
    % for each orientation Change cut
    toChangeVect = cellfun(@(x)  ~isempty(x),orientChangePtsCell);
    IDCCToChange = find(toChangeVect);  
    removeMask = zeros([ny,nx]); 
    for iChange = 1:length(IDCCToChange)
        IDC = IDCCToChange(iChange);
        % get the pixel indices to change
        pixIdx = CCInt.PixelIdxList{IDC};
        % get where to cut
        cutHere = orientChangePtsCell{IDCCToChange(iChange)};
        removeMask(pixIdx(cutHere+1))= 1; 
        pixIdx(cutHere+1) = [];
        
        CCInt.PixelIdxList{IDC} = pixIdx; % change the pixels
       % [y,x] = ind2sub(size(maxTh),pixIdx);
%         if ip.Results.TSOverlays == true; 
%         scatter(x,y,10,'r','filled');
%         text(x(1),y(1),num2str(IDC),'color','r');
%         end 
    end
   
    afterCut = zeros(size(maxTh));
    afterCut(vertcat(CCInt.PixelIdxList{:})) = 1;

    %% Reget ConnectedComponents now that have broken based on orientation
    clear CCInt csize
    
    CCInt = bwconncomp(afterCut);
    
    % for orientation calcs want to make sure to get rid of pieces less
    % than 3 ccs 
    csize = cellfun(@(x) size(x,1), CCInt.PixelIdxList);
    CCInt.PixelIdxList(csize<=3) =  [];
    CCInt.NumObjects = CCInt.NumObjects - sum(csize<=3);
    
    %GetInternal Endpoints
    
    embeddedRidgeCand1EPs = cellfun(@(x) getEndpoints(x,size(img),0,1),CCInt.PixelIdxList,'uniformoutput',0);
    %  % change EPs to pixIdx
    
    % make sure to clean out noise (ie those fragments without any endpoints)
    idxLogicNoEPs = cellfun(@(x) isempty(x) ,embeddedRidgeCand1EPs);
    % get rid of those without endpoints
    embeddedRidgeCand1EPs(idxLogicNoEPs) = [];
   % afterCut(vertcat(CCInt.PixelIdxList{idxLogicNoEPs})) = 0; % set these equal to zero in the original mask
    CCInt.PixelIdxList(idxLogicNoEPs) = [];
    CCInt.NumObjects = CCInt.NumObjects - sum(idxLogicNoEPs);
%% OPTIONAL TS Overlay  Cut-Curvature 
    if ip.Results.TSOverlays == true;
        
         
        [ny,nx] = size(img);
        TSFigs(countFigs).h  = setFigure(nx,ny,'on');
        TSFigs(countFigs).name = 'BeforeAfterCleanStepII_CutCurves';

        imshow(-img,[]);
        hold on 
        spy(beforeCut,'b',10);
        % overlay after the cut
        spy(afterCut,'r',10);
        finalIntCands = zeros([ny,nx]); 
        finalIntCands(vertcat((CCInt.PixelIdxList{:})))=1; 
        spy(finalIntCands,'g',10); 
       % finalIntCands
        
        % NOTE: eventually make an plot color coded by maxTh over the
        % img... 
        text(5,5, 'Before Curvature Cut', 'FontSize',10,'Color','b'); 
        text(5,15,'After Curvature Cut','FontSize',10,'Color','r');
        text(5,25,'Final Candidates','FontSize',10,'Color','g');  
        countFigs = countFigs +1;   % close fig
        
    
        TSFigs(countFigs).h  = setFigure(nx,ny,'on');
        TSFigs(countFigs).name = 'OrientationAndCutSites';
        beforeCut(beforeCut==0) = NaN; 
        orientMaskedByNMS = beforeCut.*rad2deg(maxTh+pi/2);
        %orientsMaskByNMS(orientsMaskByNMS==0) = NaN; 
        imagesc(orientMaskedByNMS); 
        
               cmap = get(gcf,'colormap');
               add = [0 0 0 ];
               cmap = [add;cmap]; 
               set(gcf,'colormap',cmap);
               % over lay removed 
               hold on 
               spy(removeMask,'w',10); 
    countFigs = countFigs+1; 
    end    
%% Calc displacement vectors of the ridges: will use for matching 
% NOTE 20150524 think have made this more local with getEndpoints can try
% to see if this helps matching. (DELETED 20150526 silly not necessarily
% sorted so that the vect is always facing the same way...dipshit)
     vectInt =  cellfun(@(x) [x(1,1)-x(2,1), x(1,2) - x(2,2)], embeddedRidgeCand1EPs ,'uniformoutput',0);
     dInt  = cellfun(@(x) sqrt((x(1,1)-x(2,1))^2 + (x(1,2)-x(2,2))^2),embeddedRidgeCand1EPs,'uniformoutput',0);
%     
     vectSeed =  cellfun(@(x) [x(1,1)-x(2,1), x(1,2) - x(2,2)], seedFilo1EPs ,'uniformoutput',0);
     dSeed = cellfun(@(x) sqrt((x(1,1)-x(2,1))^2 + (x(1,2)-x(2,2))^2),seedFilo1EPs,'uniformoutput',0);
%%   
    
   
    
    
    
    % 
    [intEPIdx] = cellfun(@(x) sub2ind(size(img),x(:,2),x(:,1)),embeddedRidgeCand1EPs,'uniformoutput',0);
    %
    [seedEPIdx] = cellfun(@(x) sub2ind(size(img),x(:,2),x(:,1)),seedFilo1EPs,'uniformoutput',0);
    
    % Don't need to put all pixels through matching only the ones nearest to
    % the cell body should be considered (other option would be to set up to
    % have the two end points compete.  If have time could test processing
    % time for both and see if save anything
    
    % get distTrans relative to veilStem
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
        embeddedRidgeCand1EPs(parLIdx') = [];
        
        
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
            % % get the embeddedRidgeCandCoordsClosest to the edge
            for i = 1:length(idxKeepInt)
                embeddedRidgeCand1EPsFinal{i} = embeddedRidgeCand1EPs{i}(idxKeepInt(i),:); % just take the point closest to the neurite edge
            end
            %
            %cellfun(@(x) (x(1,1)-x(1,2))^2
            embeddedRidgeCand1EPsFinal = vertcat(embeddedRidgeCand1EPsFinal{:}); % only the coords closest to the cell edge will
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
            
            embeddedRidgeCandSpur =double(labelMatCandInt1>0);
            
            % sanity
            figure
            imshow(-img,[]) 
            hold on 
            spy(filoExtSeedForInt,'r'); 
            scatter(seedFilo1EPsFinal(:,1),seedFilo1EPsFinal(:,2),10,'r','filled'); 
            quiver(seedFilo1EPsFinal(:,1),seedFilo1EPsFinal(:,2),seedFilo1EPsFinal(:,3),seedFilo1EPsFinal(:,4),0.2,'color','r')
            spy(embeddedRidgeCandSpur,'b')
            scatter(embeddedRidgeCand1EPsFinal(:,1),embeddedRidgeCand1EPsFinal(:,2),10,'b','filled'); 
            quiver(embeddedRidgeCand1EPsFinal(:,1),embeddedRidgeCand1EPsFinal(:,2),...
                embeddedRidgeCand1EPsFinal(:,3),embeddedRidgeCand1EPsFinal(:,4),0.2,'color','b'); 
            
            
            
            
            
            
            
            
            
            
            reconstruct.Int.Seed{1} = filoExtSeedForInt; % internal Seed
            reconstruct.Int.Cand{1} = embeddedRidgeCandSpur; % interal Candidates
            
 %%           
      if ip.Results.TSOverlays == true 
           
        TSFigs(countFigs).h  = setFigure(nx,ny,'on');
        TSFigs(countFigs).name = 'Before Matching';
        imshow(-img,[]); 
        hold on 
        % plot seeds 
        
        
        % plot embedded cands 
        
          
        % plot endpoints   
          
           
            

                

                
%                 imshow(img,[]);
                hold on
                
                 scatter(embeddedRidgeCand1EPsFinal(:,1),embeddedRidgeCand1EPsFinal(:,2),20,'b','filled'); % the endpoint to connect
                 scatter(seedFilo1EPsFinal(:,1),seedFilo1EPsFinal(:,2),20,'r','filled');
                 fromLabels = double(labelMatCandInt1>0); 
                 spy(embeddedRidgeCandSpur,'y'); 
                 spy(fromLabels,'b');
                 spy(filoExtSeedForInt,'r')
                 countFigs = countFigs +1; 
%                 
%                 
%             end % making trouble shoot internal figures %%
      end 
% END CLEANING 
           %% Perform the linking 
            % run through
            % maskpostconnect1 should have all the CC filos after first connection
            % between in and out - NOTE you need to keep the vect as input 
            % set up the endpoints to choose only the closest point so 
            % it makes sense not to calculate again. 
            [maskPostConnect1,linkMask1,status,TSFigs2]  = gcaConnectEmbeddedRidgeCandidatesFix(embeddedRidgeCand1EPsFinal,seedFilo1EPsFinal,filoExtSeedForInt, ...
                labelMatCandInt1,img,edgeMask,p);
            % need to get the internal Filoconnect
            
            
            
          
            reconstruct.Int.links{1} = linkMask1;
            reconstruct.Int.Seed{2} = maskPostConnect1;
            
           
        else
            maskPostConnect1 = filoExtSeedForInt;
        end % if idxKeepInt % 20140301 SEE IF YOU CAN CLEAN THIS UP
        
    else
        maskPostConnect1 = filoExtSeedForInt;
        
    end  %isempty(idxKeepInt)
    


end

