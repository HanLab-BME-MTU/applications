function [ output_args ] = GCAVeilStemTemporalImprovements(analInfo,MD,plots)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 

 
        

% check for the field to flag outliers 
if ~isfield(analInfo,'flagOutlier')
    error('Please run outlier detection');
end 

% Get all outlier frames with extra potential body pieces
outlier = arrayfun(@(x) ~isempty(x.flagOutlier),analInfo);
idxFloater = arrayfun(@(x) ~isempty(x.bodyEst.rmHiIntPieces),analInfo);
frames2Fix = find(outlier & idxFloater); 

if isempty(frames2Fix) 
    display('No Likely Truncations: No Changes Were Made'); 
else 
    % take out first and last frames 
    frames2Fix= frames2Fix(frames2Fix ~= 1 | frames2Fix ~=length(analInfo)) ;
    
for iFrame = 1:length(frames2Fix)
    
   frameC = frames2Fix(iFrame);

% load necessary pieces
floatingPieceC = analInfo(frameC).bodyEst.rmHiIntPieces;
veilStemC = analInfo(frameC).masks.neuriteEdge; 
idxEnterNeurite = analInfo(frameC).idxEnterNeurite;
[ny,nx] = size(veilStemC);

backboneMinus = analInfo(frameC-1).bodyEst.backbone;

%backbonePlus = analInfo(frameC+1).bodyEst.backbone;

%backbone = (backboneMinus | backbonePlus); 
backbone = bwmorph(backboneMinus,'thin'); 

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
 
 CCFullMask = bwconncomp(fullMask); 
 % decide later if really want to put plots here. 
 if plots == 1 
       setFigure(nx,ny,'on')
       
      imgPath = MD.getChannelPaths{1}; 
      imgPath = [imgPath filesep  MD.getImageFileNames{1}{frameC}]; % REMEMBER THEY SET THIS UP SO IN CELLS
      img  = double(imread(imgPath)); 
      
      
       imshow(-img,[])
       hold on 
       
       roiYXMaskOld = bwboundaries(veilStemC); 
       roiYXMaskNew  = bwboundaries(fullMask);
      
     
       cellfun(@(x) plot(x(:,2),x(:,1),'g','Linewidth',2),roiYXMaskNew);    
        cellfun(@(x) plot(x(:,2),x(:,1),'b','Linewidth',2),roiYXMaskOld);
       
       
       saveas(gcf,['FixFrame' num2str(frameC,'%03d') '.fig'])
       close gcf
       
end 
 
 
 
   
 if CCFullMask.NumObjects == 1
     
     analInfo(frameC).masks.neuriteEdge = fullMask; 
     analInfo(frameC).fixVeilStem = 1;
   
 else 
     analInfo(frameC).fixVeilStem =0 ; 
 end 
     
       
 


end
 
 save('analInfoTruncFix.mat','analInfo'); 
 

end 
