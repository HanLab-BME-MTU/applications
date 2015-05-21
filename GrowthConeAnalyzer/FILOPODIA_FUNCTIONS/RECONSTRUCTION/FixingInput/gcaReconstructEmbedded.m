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
    CCInt.PixelIdxList(csize<2)= [];
    CCInt.NumObjects = CCInt.NumObjects - sum(csize<2);
      
%% Clean the Seed 
    
    % keep the endpoints in a cell array that corresponds to EACH CC for Seed
    CCFiloExtSeedForInt = bwconncomp(filoExtSeedForInt);
    
    seedFilo1EPs = cellfun(@(x) getEndpoints(x,size(img)),CCFiloExtSeedForInt.PixelIdxList,'uniformoutput',0);
    
    % filter out those CCs with no or more than 2 end points
    weirdSeed = cellfun(@(x) size(x,1)~=2,seedFilo1EPs);
    seedFilo1EPs = seedFilo1EPs(~weirdSeed);
    CCFiloExtSeedForInt.PixelIdxList(weirdSeed) = [];
    CCFiloExtSeedForInt.NumObjects = CCFiloExtSeedForInt.NumObjects- sum(weirdSeed);
    
 %%  Cut Embedded Candidates At Points of High Curvature  
 % This step is required because 
 % there are very often embedded ridge candidates that are true candidates 
 % but merge with the ridge response arising at the border of the veil   
 % This helps eliminate this problem. 
   
 % Currently simply get orientations of internal ridge candidates per pixel from the
 % steerable filter output. (NOTE MB: investigate better ways to do this
 % before release)
    orient = cellfun(@(x) rad2deg(maxTh(x)+pi/2),CCInt.PixelIdxList,'uniformoutput',0);
    
    % calc gradient of orientation
    diffOrient = cellfun(@(x) abs(diff(x)),orient,'uniformoutput',0);
    
    % where orientation differences greater than 90 degrees cut segment.
    orientChangePtsCell= cellfun(@(x) find(x>20 & x <170,1,'first'),diffOrient,'uniformoutput',0);
    
    % for each orientation Change cut
    toChangeVect = cellfun(@(x)  ~isempty(x),orientChangePtsCell);
    IDCCToChange = find(toChangeVect);

    % SANITY CHECK
    beforeCut = zeros(size(maxTh));
    beforeCut(vertcat(CCInt.PixelIdxList{:})) = 1;
%%  OPTIONAL TS OVERLAY 
    if ip.Results.TSOverlays == true;
        [ny,nx] = size(img);
        TSFigs(countFigs).h  = setFigure(nx,ny,'on');
        TSFigs(countFigs).name = 'Cut_Curvature';

        imshow(-img,[]) ;
        hold on
        spy(beforeCut,'b');
        hold on
    end % figure open 
%%    
    
    for iChange = 1:length(IDCCToChange)
        IDC = IDCCToChange(iChange);
        % get the pixel indices to change
        pixIdx = CCInt.PixelIdxList{IDC};
        % get where to cut
        cutHere = orientChangePtsCell{IDCCToChange(iChange)};
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
    %
%% OPTIONAL TS Overlay  Cut-Curvature 
    if ip.Results.TSOverlays == true;
        % overlay after the cut
        spy(afterCut,'r'); 
        % NOTE: eventually make an plot color coded by maxTh over the
        % img... 
        text(5,5,'After Curvature Cut','FontSize',10,'Color','r'); 
        countFigs = countFigs +1;   % close fig
    end 
    %% Reget ConnectedComponents now that have broken based on orientation
    clear CCInt csize
    
    CCInt = bwconncomp(afterCut);
    
    % get rid of singletons
    csize = cellfun(@(x) size(x,1), CCInt.PixelIdxList);
    CCInt.PixelIdxList(csize<=2) =  [];
    CCInt.NumObjects = CCInt.NumObjects - sum(csize<=2);
    
    %GetInternal Endpoints
    
    embeddedRidgeCand1EPs = cellfun(@(x) getEndpoints(x,size(img)),CCInt.PixelIdxList,'uniformoutput',0);
    %  % change EPs to pixIdx
    
    % make sure to clean out noise (ie those fragments without any endpoints)
    idxLogicNoEPs = cellfun(@(x) isempty(x) ,embeddedRidgeCand1EPs);
    % get rid of those without endpoints
    embeddedRidgeCand1EPs(idxLogicNoEPs) = [];
   % afterCut(vertcat(CCInt.PixelIdxList{idxLogicNoEPs})) = 0; % set these equal to zero in the original mask
    CCInt.PixelIdxList(idxLogicNoEPs) = [];
    CCInt.NumObjects = CCInt.NumObjects - sum(idxLogicNoEPs);

%% Calc displacement vectors of the ridges: will use for matching 
    vectInt =  cellfun(@(x) [x(1,1)-x(2,1), x(1,2) - x(2,2)], embeddedRidgeCand1EPs ,'uniformoutput',0);
    dInt  = cellfun(@(x) sqrt((x(1,1)-x(2,1))^2 + (x(1,2)-x(2,2))^2),embeddedRidgeCand1EPs,'uniformoutput',0);
    
    vectSeed =  cellfun(@(x) [x(1,1)-x(2,1), x(1,2) - x(2,2)], seedFilo1EPs ,'uniformoutput',0);
    dSeed = cellfun(@(x) sqrt((x(1,1)-x(2,1))^2 + (x(1,2)-x(2,2))^2),seedFilo1EPs,'uniformoutput',0);
%%   
    
   
    
    
    
    
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
            reconstruct.Int.Seed{1} = filoExtSeedForInt; % internal candidates
            reconstruct.Int.Cand{1} = embeddedRidgeCandSpur; % interal Candidates
            
            
            
            
%% Troubleshoot figure internal
%             
%             if ip.Results.TSOverlays == true
%                 TSFigs(countFigs).h = setFigure(nx,ny,'on'); 
%                 TSFigs(countFigs).name = 'beforeAfterSpur';
%                
%                 
%                 imshow(img,[]);
%                 hold on
%                 
%                 spy(embeddedRidgeCand,'g'); % green is the original filo detection after thresholding 
%                 spy(embeddedRidgeCandSpur,'b'); % blue is after spur
                
               
             
%                 countFigs = countFigs +1; 
                

                
%                 imshow(img,[]);
%                 hold on
                
%                 scatter(embeddedRidgeCand1EPsFinal(:,1),embeddedRidgeCand1EPsFinal(:,2),20,'b','filled'); % the endpoint to connect
%                 scatter(seedFilo1EPsFinal(:,1),seedFilo1EPsFinal(:,2),20,'r','filled');
%                 spy(embeddedRidgeCandSpur,'b');
%                 countFigs = countFigs +1; 
%                 
%                 
%             end % making trouble shoot internal figures %%
% END CLEANING 
           %% Perform the linking 
            % run through
            % maskpostconnect1 should have all the CC filos after first connection
            % between in and out - NOTE you need to keep the vect as input 
            % set up the endpoints to choose only the closest point so 
            % it makes sense not to calculate again. 
            [maskPostConnect1,linkMask1,status]  = gcaConnectEmbeddedRidgeCandidates(embeddedRidgeCand1EPsFinal,seedFilo1EPsFinal,filoExtSeedForInt,labelMatCandInt1,vectSeed,vectInt,dSeed,dInt,p);
            % need to get the internal Filoconnect
            
            
            
            if ip.Results.detectEmbedded == true;
                if status ==1
                    spy(linkMask1,'y',10); % plot the link mask
                end
                TSFigs(countFigs).name = 'graphConnect'; 
                TSFigs(countFigs).h = graphConnectFig; 
                countFigs = countFigs +1; 
                
                %saveas(gcf,[linkDir filesep 'internalInputEPsWithLinks' num2str(iFrame) '.tif']);
            end
            reconstruct.Int.links{1} = linkMask1;
            reconstruct.Int.Seed{2} = maskPostConnect1;
            
           
        else
            maskPostConnect1 = filoExtSeedForInt;
        end % if idxKeepInt % 20140301 SEE IF YOU CAN CLEAN THIS UP
        
    else
        maskPostConnect1 = filoExtSeedForInt;
        
    end  %isempty(idxKeepInt)
    


end

