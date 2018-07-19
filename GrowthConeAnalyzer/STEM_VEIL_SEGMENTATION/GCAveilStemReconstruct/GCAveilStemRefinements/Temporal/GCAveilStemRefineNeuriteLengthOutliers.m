function [ analInfo ] = GCAVeilStemTemporalImprovements(analInfo)
%GCAVeilStemTemporalImprovments: Currently a simple function that attempts
%to avoid truncations in the segmentation by identfying outliers in the
%neurite length time series.  If this frame likewise has an unused body
%segment it is highly probable that this body part should be connected. 
% the easiest way to reconnect is to simply use the reconstructed backbone
% from the surrounding frames- as this backbone tends to change little
% within small intervals. This is kept as a second step as movies with
% low time resolution or undergoing more dramatic dynamic behavior may not
% be able to use this prior to aid in segmentation. 

% check for the field to flag outliers 
if ~isfield(analInfo,'flagOutlier')
    error('Please run outlier detection');
end 

% Get all outlier frames with extra potential body pieces
outlier = arrayfun(@(x) ~isempty(x.flagOutlier),analInfo);
idxFloater = arrayfun(@(x) ~isempty(x.rmHiIntPieces),analInfo);
frames2Fix = find(outlier & idxFloater); 

if isempty(frames2Fix) 
    display('No Likely Truncations: No Changes Were Made'); 
else 
    % take out first and last frames 
    frames2Fix= frames2Fix(frames2Fix ~= 1 | frames2Fix ~=length(analInfo)) ;
    
for iFrame = 1:length(frames2Fix)
    display(['Fixing Frame ' num2str(frames2Fix(iFrame))]); 
   frameC = frames2Fix(iFrame);

% load necessary pieces
floatingPieceC = analInfo(frameC).rmHiIntPieces;
veilStemC = analInfo(frameC).finalMask; 

idxEnterNeurite = analInfo(frameC).idxEnterNeurite;
[ny,nx] = size(veilStemC);
framesForAndBack = frameC-5:frameC+5; % need a flag for early/late  guys 
backbones = arrayfun(@(x) analInfo(x).backbone,framesForAndBack,'uniformoutput',0);
            sumBB = zeros(ny,nx);
            
            for iBB = 1:numel(backbones)
                sumBB = sumBB+backbones{iBB}; % likely better way than a loop for this but quick fix
            end
            % find pixels in structure that were more in more than 5 frames (to remove
            % outliers)
            sumBB(sumBB<=2) =0;
            sumBB(sumBB>0) = 1;% make a logical mask
            backbone = bwmorph(sumBB,'thin'); 

%backboneMinus = analInfo(frameC).bodyEst.backbone;

%backbonePlus = analInfo(frameC+1).bodyEst.backbone;

%backbone = (backboneMinus | backbonePlus); 
%backbone = bwmorph(backboneMinus,'thin'); 

putTogether = (floatingPieceC | veilStemC | backbone); 

newBodyMask = (floatingPieceC | veilStemC); 
notBody = double(backbone).*~double(newBodyMask);

% 12-11 note just make sure to thin it here. check why this thin step was
% miseed before ...though I previously implemented it... 

notBody = bwmorph(notBody,'thin','inf'); 

newBodyXY = bwboundaries(newBodyMask); 
         idxBody = cellfun(@(x) sub2ind(size(newBodyMask),x(:,1),x(:,2)),newBodyXY,'uniformoutput',0);
         idxBody = vertcat(idxBody{:}); 
         bodyNoFill = zeros([ny,nx]); 
         bodyNoFill(idxBody) = 1; 

 
 %% 2013_12_08 check this - in theory should not have junctions if did this cleanly (rethink) 
% Break junctions 
% also don't want to break junctions in the original backbone... 
  nn = padarrayXT(double(notBody~=0), [1 1]);
  sumKernel = [1 1 1];
  nn = conv2(sumKernel, sumKernel', nn, 'valid');
  forJunctBreak = notBody~=0; 
  % for now just make sure not f-ing up the original backbone 
  % 
  %forJunctBreak(backboneInfo(iFrame).backboneSeedMask==1) = 0; 
 test = notBody|bodyNoFill;% get the backbone pieces and the no fill body for erosion 
        test2 = bwmorph(test,'thin',inf); 
        backbone(test & ~test2) = 0 ; % not  sure what I was doing here 
        %get rid of singletons : this might not be necessary if don't break
        %junctions in the first step 
        CCTest2 = bwconncomp(test2); 
        csize = cellfun(@(x) length(x),CCTest2.PixelIdxList); 
        backbone(vertcat(CCTest2.PixelIdxList{csize==1})) = 0; % set these = to zero   
        CCTest2.PixelIdxList(csize==1) = []; 
        CCTest2.NumObjects = CCTest2.NumObjects - sum(csize==1);
        test2= labelmatrix(CCTest2);
     test2 = test2>0;
        [EPs,~,coords] = skel2graph2D(test2);

  if ~isempty(vertcat(coords{:}))
            indEP = sub2ind([ny,nx],EPs(:,1),EPs(:,2));
         
           
        
            
            %idxSave = find(indEP == idxEnterNeurite); %% note sometimes bug here... should make so reiterate if this fails...
            %instead of doing find use intersect
            overlap = intersect(idxEnterNeurite,indEP); 
            % save the entering neurite pieces from erosion
            if ~isempty(overlap)
                 % maybe add to take it up or down a scale...but for now
                 % just leave it
            idxSave = arrayfun(@(x) find(indEP == overlap(x)),1:length(overlap)); 
            idxSave = idxSave'; 
           
                EPs(idxSave,:) = [];
                coords(idxSave) = [];
            end
            
            
            
            if ~isempty(vertcat(coords{:}));
                coordBBOver = vertcat(coords{:});
                idxBBOver = sub2ind([ny,nx],coordBBOver(:,1),coordBBOver(:,2));
                idxEP = sub2ind([ny,nx],EPs(:,1),EPs(:,2));
                backbone([idxBBOver ; idxEP]) = 0;
            end
        end
         
 
  %%      
        backbone2Dil = backbone;
        backbone2Dil(backbone==1 & newBodyMask==1) = 0; 
        % for now the dilation might be the best 
        
  
  %%% FOR NOW JUST ARBITRARILY RE-DILATE
            dilBB = imdilate(backbone2Dil,strel('disk',4));
            %          else
  
 fullMask = dilBB | newBodyMask;
 figure; 
 imshow(fullMask,[]) ; 
 
 CCFullMask = bwconncomp(fullMask); 
   %% added 20140104
%   roiYXAll = bwboundaries(fullMask);
%     pixIndicesAll = sub2ind(size(fullMask),roiYXAll{1}(:,1),roiYXAll{1}(:,2));
%     
%     
% 
%     % load the old thin body pixels
%     pixIndicesThinBody =  analInfo(frameC).bodyEst.pixIndThinBody;
%     % anything for now that is not thin body is new thick body
%     pixIndicesThickBodyNew = setdiff(pixIndicesAll,pixIndicesThinBody);
%     
%     thickBodyMask = zeros(size(fullMask));
%     thickBodyMask(pixIndicesThickBodyNew)=1;
%     thickBodyMask = logical(thickBodyMask);
%     % record
%     analInfo(frameC).masks.thickBodyMask = thickBodyMask;
%     analInfo(frameC).bodyEst.pixIndThickBody = pixIndicesThickBodyNew;
    
%% 
 
 
 
 if CCFullMask.NumObjects == 1
     analInfo(frameC).masks.neuriteEdge = fullMask; 
     
     analInfo(frameC).outlierFixFlag = 1; 
 else 
     analInfo(frameC).outlierFixFlag =0 ; 
 end 
     
       
 


end
 
 % in the future make a new folder and save 
 save('analInfoTruncFix.mat','analInfo'); 
 
 
 

end 
