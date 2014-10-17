function [ filoInfo ] = gcaRecordFilopodiaInformation( CCFiloObjs,img,maxRes,maxTh,edgeMask,bodyMask,analInfoC,normalsC,smoothedEdgeC)
% gcaRecordFilopodiaInformation: documents information such as pixel ID from 
% base to tip, orientation, etc, about each of the filopodia objects input as connected components. 
% internal function of Growth Cone Analyzer
%
% NOTE TO ME: PLEASE FIX INPUT BEFORE RELEASE
%
% INPUT: 
%       CCFiloObjs: connected component output of all possible filopodia to
%       to be documented
%       
%       img (needed for intensity information)
%      
%       maxTh, maxRed (consolidate from analInfo) 
%       edgeMask,bodyMask (same thing)  you can potentially consolidate
%       
%       analInfoC: structure containing the current analInfo for iFrame,
%       here is where you have most of your information regarding the
%       previous imageProcessing
%
%       normalsC: the normal corresponding to iFrame, output of sam's 
%       protrusion bundle. if empty don't calculate filopodia orientations
%
%       smoothedEdgeC: the spline fit edge estimation from Sam's protrusion
%       software. if empty don't calculate filopodia orientations
%       
% OUTPUT: 
%      filoInfo: an 1 x nFilo structure with separate fields containing all relavent information
%      for each filopodia in current frame. Designed to be somewhat
%      adaptable (more or less information can potentially be stored here)
%      Fields are designed here to be somewhat descriptive. 
%      Internal and External Filo information are stored in separate fields
%      starting with Int_ or Ext_ respectively      
%

countFilo = 1; 
[ny,nx] = size(maxTh); 
dims = [ny,nx]; 
filoInfo = struct([]); 
% set up filopodia structure from the get go so don't have trouble with
% structure having different number of fields if don't fit internal
% filopodia (therefore default will be empty) 
toAdd{1} = 'Ext_'; 
toAdd{2} = 'Int_'; 

fields{1} = 'coordsXY'; 
fields{2} = 'dirAtTip'; 
fields{3} = 'endpointCoord'; 
fields{4} = 'intensities'; 
fields{5} = 'pixIndices'; 
fields{6} = 'pixIndicesBack'; 
fields{7} = 'pixIndicesFor'; 
fields{8} = 'response'; 
fields{9} = 'vectFilo'; 
fields{10} = 'maskIndices'; 

for i = 1:2
    for j = 1:numel(fields)
        filoInfo(1).([toAdd{i} fields{j}]) = [];
    end
end
filoInfo.cross = []; 
filoInfo.type = []; 
filoInfo.groupCount = []; 
filoInfo.bodyType = []; 
%if ~isempty(normalsC) % always initiate
filoInfo.orientation = []; 
filoInfo.localVectAttach = []; 
filoInfo.localVectFilo = []; 
%end 
filoInfo = orderfields(filoInfo); 

for iFiloObj = 1:numel(CCFiloObjs.PixelIdxList)
    % make the individual mask for each filopodia simplies the labeling and
    % makes more intuitive: if time there might be a more clever way to save
    % time on this computationally and not make it so much of a loop. But I
    % found it helps at least now for troubleshooting and information
    % organization
    maskCurrent = zeros(dims);
    maskCurrent(CCFiloObjs.PixelIdxList{iFiloObj})=1;
    
    
    maskCurrentInt = maskCurrent.*bodyMask;
    maskCurrentExt = maskCurrent.*~bodyMask;
    if sum(maskCurrentExt(:))>0 && sum(maskCurrentInt(:))>0; % if external filo check for internal filo
        intFlag =2; % fit response have
    elseif sum(maskCurrentExt(:))>0 && sum(maskCurrentInt(:))==0;
        intFlag =1; % just project forward to search for missing response % this option seems to be
        % VERY noisymatlab
    else
        intFlag =0 ;
    end
    
    maskCurrentInt = maskCurrentInt|edgeMask;
    maskCurrentExt = maskCurrentExt | edgeMask;
    test = 1; 
    if test == 1
        imshow(img,[]) 
        hold on 
        spy(maskCurrentInt,'b'); 
        spy(maskCurrentExt,'r'); 
       
    end 
%%%% GET EXTERNAL INFORMATION FIRST %%%% 
    
       [verticesEP,verticesBP, edgePathCoord] = skel2graph2D(maskCurrentExt);
       
       if ~isempty(verticesEP)
       
%        
           x = walkFiloForAndBack([],verticesEP,edgePathCoord,maxTh,maxRes,img,0,10); % typiclly use 10 pixels forward for the external 
           
           
           if length(verticesEP(:,1)) > 1
               % means that there is something fishy going on in theory
               % here set it up to only have 1 endpoint...
               % this can happen that the external filopodia are connected
               % by a strong internal response along edge...
               % and need to potentially cut this also is potentially an
               % indication that you have bad body segmentation % Consider this
               % a QUICK FIX for now
               % for now just loop through the x info
               for i = 1:length(verticesEP) 
                   x(i).cross = 0;
                   x(i).type =0;
                   x(i).groupCount = countFilo; 
                   x(i).bodyType = NaN;
                   if ~isempty(normalsC); 
                   x(i).orientation = NaN; % for now just leave blank might want to do separately
                   x(i).localVectFilo = NaN;
                   end
                   fieldsx = fieldnames(x(i)); 
                   if exist('filoInfo','var')
                       fieldsFilo = fieldnames(filoInfo);
                       fieldsToAdd = setdiff(fieldsFilo,fieldsx);
                       for iField = 1:length(fieldsToAdd)
                           x(i).(char(fieldsToAdd(iField))) = NaN; % set all internal parameters to NaN; can maybe remove now
                       end
                     
                       filoInfo = orderfields(filoInfo);
                   end
                   
                   toSave = orderfields(x(i));
                   
                  filoInfo(countFilo) = toSave;
                   countFilo = countFilo+1; 
                   if exist('toSave','var')
                   clear toSave
                   end
               end
               clear x 
               intFlag = 0; % don't continue the internal likely a mess.
               
               
               
           else % if normal set these fields 
               
               
     
               x.cross = 0; % do not mark as cross
               x.type = 0; % type is such that it is attached to the body 
               x.groupCount = countFilo; % this will initiate the branch counting
           pathCoords = edgePathCoord{:}; %    
             
          baseFilo = [pathCoords(end,2),pathCoords(end,1)];
%% GET LOCAL INFORMATION SURROUNDING FILOPODIA %%%%
          %% Thick body thin body test 
           
          % branchpoints = bwmorph(maskCurrentExt,'branchpoints');
          % make mask of the end of path and dilate to get the surrounding
          % area 
          testMask = zeros(size(img)); 
       
          testMask(pathCoords(end,1),pathCoords(end,2)) = 1; 
          % find the neighborhood 
          testMask = imdilate(testMask,strel('disk',2));
           idx = find(testMask==1) ;
          
          % NOTE: Maria you had a bug here in your neuriteEstimate software
          % when fixing your dilation problem in Control 06 Test SetII 
          % where you forgot to save the pixIndThickBody and pixIndThinBody
          % it is easy however to calculate so just fix here. 
          if ~isfield(analInfoC.bodyEst,'pixIndThickBody'); 
              % add 
              thickBodyMask = logical(analInfoC.masks.thickBodyMask); 
              analInfoC.bodyEst.pixIndThickBody = find(thickBodyMask==1); 
              neuriteEdgeMask = analInfoC.masks.neuriteEdge; 
              thinBodyMask = neuriteEdgeMask.*~thickBodyMask;  
              analInfoC.bodyEst.pixIndThinBody = find(thinBodyMask==1); 
          end 
          
           
           test1 = intersect(idx,analInfoC.bodyEst.pixIndThickBody);
           test2 = intersect(idx,analInfoC.bodyEst.pixIndThinBody);
           
            thickBody = ~isempty(test1);
           thinBody = ~isempty(test2);
           if thickBody + thinBody ==2
               lengths = [length(test1) length(test2)];
               tiebreaker = find(lengths==max(lengths));
               if tiebreaker ==1
                   thinBody = 0;
               else
                   thickBody = 0;
               end
           end
           
           % set body type
           if thickBody == 1 && thinBody==0
               bodyType = 1; % thick is 1 ;
           elseif thickBody  == 0 && thinBody ==1
               bodyType  = 2;
           else % just in case something weird
               bodyType = NaN;
           end
            x.bodyType = bodyType;
                clear testMask    
%                      
%% ORIENTATION MEASUREMENT %%  
           % Get local normal vectors from protrusion output
           % (KD tree here is easiest as it maintains indexing of
           % protrusion output- NOTE: might want to change above thick/thin
           % region identification.
           if ~isempty(normalsC) % if have run through the protrusion software
               % get edge coordinates within 3 pixels of the base of the
               % filo using the KD tree. 
               
               [idx, dist] = KDTreeBallQuery(smoothedEdgeC,baseFilo,3); 
               idx = idx{:};
            
                   
               % Note can either get orientation of filo from the maxTh output of the steerable filter or from just calculating
               % a small local vector. we will see which one is cleaner ... so far I tend
               % to favor the small vector...
               avgNormLocal = mean(normalsC(idx,:));% might want to change to a majority? (i don't think these vectors are really normalized)
               % test length of pathCoords
               pixFilo = size(pathCoords,1);
               if pixFilo <= 4
                   back = pixFilo-1;
               else
                   back = 4;
               end
               % Local vector calc (NOTE could also potentially get this angle from the
               % maxTh data usually take 3-4 pixels "back" should really call for...
               % pathCoords(end) is where the filo intersects with body so vector is in
               % the direction away from the cell body
               if back ~=0  % if the filo is long enough and you are averaging over enough of the edge proceed 
                   localVectFilo = [pathCoords(end-back,2)-pathCoords(end-1,2),pathCoords(end-back,1)-pathCoords(end-1,1)];
                   vectLength = sqrt((pathCoords(end-back,2)-pathCoords(end-1,2)) ^2 + (pathCoords(end-back,1) - pathCoords(end-1,1))^2);
               else 
                   localVectFilo = [NaN,NaN];
                   vectLength = NaN; 
               
               end % back ~=0 
               if length(idx) >=3 % if there is not enough local edge pixels over which to average (this is very common when put only part of the body in the protrusion software)
               
                   normLength = sqrt(avgNormLocal(1)^2 + avgNormLocal(2)^2);% check this... 08-10
               else 
                   normLength = NaN; 
                   avgNormLocal = [NaN NaN]; 
               end % idx > = 3 
               
               % calculate angle to body 
                   cosAngle = dot(avgNormLocal,localVectFilo)/vectLength/normLength;
                   angle = rad2deg(acos(cosAngle));
                   angleToBody = 180- angle -90;
                 
               x.orientation = angleToBody;
               x.localVectAttach = avgNormLocal; % for now just save the normal vector 
               x.localVectFilo= localVectFilo;
               else 
                   x.orientation = []; % keep empty to maybe calculate later... ( the protrusion data was not run) 
                   x.localVectAttach = []; 
                   x.localVectFilo = []; 
               
           end % isempty normalC    
           end 
       else intFlag =0; % don't countinue if verticesEP is empty
       end % verticesEP
    clear verticesEP edgePathCoord
    %% 2013_07_14 note think this was the old way of doing things before internal 
    % was done via graph match..should likely take out this option ..all
    % internal should have a corresponding external by the way they were
    % saved. ...
    
    % OLD NOTES : for now if there is no corresponding external filopodia there is a high chance
    % that  the internal signal is just noise we will not record % 
    
    if intFlag ==2 % fit internal filo using response % 2013_07_14 again should take out before release....
        % get internal info
        [verticesEP,~, edgePathCoord] = skel2graph2D(maskCurrentInt);
       
        if length(verticesEP(:,1))==1
            x = walkFiloForAndBack(x,verticesEP,edgePathCoord,maxTh,maxRes,img,1,20,bodyMask );% 2013_07_14 try 15 pixels 
            % or if just start fitting to noise. 
           
                x = orderfields(x);
                filoInfo(countFilo) = x;
                clear x
                countFilo = countFilo +1;
        else % don't record anything likely wack
            % always safer to record the ext but not the int  if suspicious
            % 
            % 
            %intFlag =1; 
            %ADDED 20141017
              fieldsx= fieldnames(x); 
                       fieldsFilo = fieldnames(filoInfo);
                       fieldsToAdd = setdiff(fieldsFilo,fieldsx);
                       for iField = 1:length(fieldsToAdd)
                           x.(char(fieldsToAdd(iField))) = NaN; % set all internal parameters to NaN; can maybe remove now
                       end
                      
                       
          filoInfo = orderfields(filoInfo);          
          toSave = orderfields(x); 
        filoInfo(countFilo) = toSave;
        
        clear x
        countFilo = countFilo+1;
            
            
            
            
            
          
        end
%         if isempty(verticesEP(:,1)) % check to see if EP empty: if still  switch to internal flag 
%             
%             intFlag =1; % something fishy is up with internal so just ignore but make sure to include the external 
%         end 
%             
       
        
    end % if intFlag ==2 
    if intFlag ==1  % try to fit it back anyway and see what get this is typically very noisy
%         coordsXY = x.Ext_coordsXY(1:2,:);
%         verticesEP = [coordsXY(1,2) coordsXY(1,1)]; % this input yx
%         
%         edgePathCoord{1} = [[coordsXY(2,2), coordsXY(2,1)];[coordsXY(1,2), coordsXY(1,1)]];
%         x = walkFiloForAndBack(x,verticesEP,edgePathCoord,maxTh,maxRes,img,1);
%        x = orderfields(x); 
                       
                       fieldsx= fieldnames(x); 
                       fieldsFilo = fieldnames(filoInfo);
                       fieldsToAdd = setdiff(fieldsFilo,fieldsx);
                       for iField = 1:length(fieldsToAdd)
                           x.(char(fieldsToAdd(iField))) = NaN; % set all internal parameters to NaN; can maybe remove now
                       end
                      
                       
          filoInfo = orderfields(filoInfo);          
          toSave = orderfields(x); 
        filoInfo(countFilo) = toSave;
        
        clear x
        countFilo = countFilo+1;
        
    else % dont' count filo and don't record likely noise
        
        
    end
end



% totalMaskInt = zeros(dims); 
% forMaskIn = vertcat(filoInfo(:).Int_pixIndices); 
% forMaskIn = forMaskIn(~isnan(forMaskIn)); 
% totalMaskInt(forMaskIn)=1; % make sure to take out any NaN here. 
% totalMaskInt = totalMaskInt |edgeMask; 
  filoInfo = fitLinescansNew(filoInfo,[ny,nx],0,11,1,0); % for now don't plot
  
%   totalMaskExt = zeros(dims); 
%   totalMaskExt(vertcat(filoInfo(:).Ext_pixIndices))=1; 
%   totalMaskExt = totalMaskExt|edgeMask; 
  filoInfo = fitLinescansNew(filoInfo,[ny,nx],0,11,0,0); %for now don't plot
  
  % calculate filopodia info into out;
  for iFilo = 1: numel(filoInfo)
      
      if ~isnan(filoInfo(iFilo).Int_length) 
      filoInfo(iFilo).totalLength = filoInfo(iFilo).Int_length + filoInfo(iFilo).Ext_length; 
      filoInfo(iFilo).percentInternal = filoInfo(iFilo).Int_length./filoInfo(iFilo).totalLength; 
      else 
          filoInfo(iFilo).totalLength = filoInfo(iFilo).Ext_length;
          if ~isnan(filoInfo(iFilo).Ext_length)
          filoInfo(iFilo).percentInternal = 0; 
          else 
              filoInfo(iFilo).percentInternal = NaN; 
          end 
      end 
      
  end 
  
  
  
