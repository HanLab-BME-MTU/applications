function [ filoInfo ] = gcaRecordFilopodiaSeedInformation( CCFiloObjs,img,filoBranchC,edgeMask,veilStemMaskC,varargin)
% gcaRecordSeedFilopodiaInformation: documents information about strong
% filopodia ridge candidates intially connected to the veil stem estimation.
% (the seed for the next reconstruction step)
% If the embedded reconstruct is run previously: it will likewise
% document the information corresponding to these internal ridges that
% were connected in this step and that havehigh probability of connecting to the high confidence
% external filopodia.
% Information currently includes the ordered path of the filopodia pixel
% coordinates (to be used for fitting) and the orientation relative to the
% veil/stem about each of the filopodia objects (if protrusion vector
% calculations have been run).
% Note currenty the internal filopodia reconstruction only uses seeds from
% the highest confidence filopodia ie those connected directly to the veil
% stem estimate- the reasoning is simply that the ridge response at the veil/stem junction
% is likely relatively strong for embedded filopodia and therefore
% lower confidence filopodia attached in later iterations are much less
% likely to have a viable internal counterpart.
%
% Internal function of Growth Cone Analyzer
%
%
%% INPUT
%
%       CCFiloObjs: connected component output of all possible filopodia to
%       to be documented
%
%       img (needed for intensity information)
%
%     
%       edgeMask,veilStemC (same thing)  you can potentially consolidate
%
%       filoBranchC: structure with fields.
%       .filterInfo (output from small scale steerable filter)
%           .maxTh: rxc double of orientation information
%           .maxRes: rxc double of response (pre-NonMaximumSuppression)
%           .scaleMap: rxc double of scale
%           .ThreshNMS: rxc double NMS response after thresholding
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
%% INPUT
%% Check Input
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('CCFiloObjs');
ip.addRequired('img');
ip.addRequired('filoBranchC',@isstruct);
ip.addRequired('edgeMask'); 
ip.addRequired('veilStemMaskC'); 

% Optional For Normal Calculations
ip.addOptional('normalsC',[]);
ip.addOptional('smoothedEdgeC',[]); 

ip.addParameter('numPixSearchForward',10,@isscalar)
ip.addParameter('numPixSearchForwardEmbed',20,@isscalar);  


%   
% ip.addParameter('maxRadiusLinkFiloOutsideVeil',10);
% ip.addParameter('minCCRidgeOutsideVeil',3);


ip.parse(CCFiloObjs,img,filoBranchC,edgeMask,veilStemMaskC,varargin{:});
p = ip.Results;
p = rmfield(p,{'CCFiloObjs','img','filoBranchC','edgeMask','veilStemMaskC'});







%% Initiate 
countFilo = 1;

[ny,nx] = size(img);
dims = [ny,nx];

maxTh = filoBranchC.filterInfo.maxTh;
maxRes = filoBranchC.filterInfo.maxRes; 

normalsC = ip.Results.normalsC; 
smoothedEdgeC = ip.Results.smoothedEdgeC; 

% Initiate output
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
%filoInfo.bodyType = [];
%if ~isempty(normalsC) % always initiate
filoInfo.orientation = [];
filoInfo.localVectAttach = [];
filoInfo.localVectFilo = [];
%end
filoInfo = orderfields(filoInfo);

%% Start recording information 
for iFiloObj = 1:numel(CCFiloObjs.PixelIdxList)
    % make the individual mask for each filopodia simplies the labeling and
    % makes more intuitive: if time there might be a more clever way to save
    % time on this computationally and not make it so much of a loop. But I
    % found it helps at least now for troubleshooting and information
    % organization
    
    maskCurrent = zeros(dims);
    maskCurrent(CCFiloObjs.PixelIdxList{iFiloObj})=1;
    
    
    % Check for For Embedded Actin Content
    maskCurrentInt = maskCurrent.*veilStemMaskC;
    maskCurrentExt = maskCurrent.*~veilStemMaskC;
    
    if sum(maskCurrentExt(:))>0 && sum(maskCurrentInt(:))>0; % if external filo check for internal filo
        intFlag =2; % Embedded Found Fit Response 
    elseif sum(maskCurrentExt(:))>0 && sum(maskCurrentInt(:))==0;
        intFlag =1; % Filo Detected Outside Veil Only: No Embedded Fitting Needed 
    else
        intFlag =0 ; % Something Fishy...Don't record anything, Could maybe happen if either piece is very small
        % primarily also part of the veil border 
    end
    
    maskCurrentInt = maskCurrentInt|edgeMask;
    maskCurrentExt = maskCurrentExt | edgeMask;
%     test = 1;
%     if test == 1
%         imshow(img,[])
%         hold on
%         spy(maskCurrentInt,'b');
%         spy(maskCurrentExt,'r');
%         
%     end
%% GET EXTERNAL INFORMATION FIRST %%%%
    
    [verticesEP,verticesBP, edgePathCoord] = skel2graph2D(maskCurrentExt);
    
    if ~isempty(verticesEP)
        %%
        % To Replace. add the parameter 
        % x = gcaProjectAndRecordFiloCoords(verticesEP,edgePathCoord,maxTh,maxRes,img,'ip.Results.numPixSearchForward',10);
        %
        
        x = walkFiloForAndBack([],verticesEP,edgePathCoord,maxTh,maxRes,img,0,10); % typiclly use 10 pixels forward for the external
        %%
        
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
                %x(i).bodyType = NaN;
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
%             %% Thick body thin body test
            
            % branchpoints = bwmorph(maskCurrentExt,'branchpoints');
            % make mask of the end of path and dilate to get the surrounding
            % area
            testMask = zeros(size(img));
            
            testMask(pathCoords(end,1),pathCoords(end,2)) = 1;
            % find the neighborhood
            testMask = imdilate(testMask,strel('disk',2));
            idx = find(testMask==1) ;
 %% Take out the thick/thin body designation           
            % NOTE: Maria you had a bug here in your neuriteEstimate software
            % when fixing your dilation problem in Control 06 Test SetII
            % where you forgot to save the pixIndThickBody and pixIndThinBody
            % it is easy however to calculate so just fix here.
            %if ~isfield(veilStemMaskC,'pixIndThickBody');
                % add
%                 thickBodyMask = logical(analInfoC.masks.thickBodyMask);
%                 analInfoC.bodyEst.pixIndThickBody = find(thickBodyMask==1);
%                 neuriteEdgeMask = analInfoC.masks.neuriteEdge;
%                 thinBodyMask = neuriteEdgeMask.*~thickBodyMask;
%                 analInfoC.bodyEst.pixIndThinBody = find(thinBodyMask==1);
%             end
%             
            
%             test1 = intersect(idx,analInfoC.bodyEst.pixIndThickBody);
%             test2 = intersect(idx,analInfoC.bodyEst.pixIndThinBody);
%             
%             thickBody = ~isempty(test1);
%             thinBody = ~isempty(test2);
%             if thickBody + thinBody ==2
%                 lengths = [length(test1) length(test2)];
%                 tiebreaker = find(lengths==max(lengths));
%                 if tiebreaker ==1
%                     thinBody = 0;
%                 else
%                     thickBody = 0;
%                 end
%             end
%             
% %             % set body type
%             if thickBody == 1 && thinBody==0
%                 bodyType = 1; % thick is 1 ;
%             elseif thickBody  == 0 && thinBody ==1
%                 bodyType  = 2;
%             else % just in case something weird
%                 bodyType = NaN;
%             end
%             x.bodyType = bodyType;
%             clear testMask
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
 %% SANITY CHECK 
                sanityCheck =0; 
                if sanityCheck == 1
                if iFiloObj == 1
                    
                   imshow(-img,[]); 
                   hold on 
                end
                   roiYX = bwboundaries(veilStemMaskC); 
                   cellfun(@(x) plot(x(:,2),x(:,1),'b'),roiYX); 
                   scatter(smoothedEdgeC(idx,1),smoothedEdgeC(idx,2),10,'b','filled'); 
                   quiver(smoothedEdgeC(idx(3),1),smoothedEdgeC(idx(3),2), avgNormLocal(1),avgNormLocal(2),10,'filled','color','c','Linewidth',2); 
                   quiver(smoothedEdgeC(idx,1),smoothedEdgeC(idx,2),normalsC(idx,1),normalsC(idx,2),'color','g'); 
                end
%%                
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
                    %% 
                    if sanityCheck == 1
                        scatter(pathCoords(end-back:end,2),pathCoords(end-back:end,1),10,'filled','r');
                        quiver(pathCoords(end,2),pathCoords(end,1),localVectFilo(1),localVectFilo(2),10,'filled','r'); 
                    end
                    %%
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
                cosAngle = dot(avgNormLocal(1:2),localVectFilo)/vectLength/normLength;
                angleToBody = acosd(cosAngle);
                %%
                if sanityCheck ==1 
                    text(pathCoords(end,2),pathCoords(end,1),num2str(angleToBody,3),'color','k'); 
%                     saveas(gcf,[num2str(iFiloObj,'%03d') '.png']); 
%                    close gcf 
                end
                %%
                % angleToBody = 180- angle -90;
                
                x.orientation = angleToBody; % in degrees
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
            %% To Replace 
            x = walkFiloForAndBack(x,verticesEP,edgePathCoord,maxTh,maxRes,img,1,20,veilStemMaskC );% 2013_07_14 try 15 pixels
            % or if just start fitting to noise.
            %x =
            %gcaProjectAndRecordFiloCoords(verticesEP,edgePatchCoord,maxTh,maxRes,img,veilStemMaskC, ...
            %'embeddedFlag',true,'numPixSearchForwardEmbed',p.numPixSearchForwardEmbed);
            %
            %%
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
        
        
        
    end % if intFlag ==2
    if intFlag ==1  % try to fit it back anyway and see what get this is typically very noisy
        
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
    count = 0; 
    if intFlag~=2; 
        display(num2str(intFlag)); 
        count = count+1;
    end 
     
end
if sanityCheck==1 
    saveas(gcf,'test.fig'); 
    saveas(gcf,'test.png'); 
    close gcf
end 





