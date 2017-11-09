function [ subRois,xVect,yVect ] = GCAgetGrowthConeSubRegions(frontPixelInfo,angle,imageSize,vectLength,roiMask)
%GCAgetGrowthConeSubRegions: This function makes a front subregion to the growth cone 
% finding the largest distance transform along a line of pixels from the 
% veil/stem estimated (calculated in previous step), finding the local 
% direction of the skeleton in this region, and adding +/- 45 directional 
% degrees to this region in order to automate a subroi. 
% 
% frontPixels: n x 2 double array of pixel n = the number of pixels as defined by 
%              the distCutOff in calculateDistance
%              column 1 is the pixel indices and column 2 is the distance 
%              transform from the veil/stem estimation (large protrusions)
%     
% angle: scalar defining the angle (in degrees( to either side of the front axis is
%        included in the front regional parameter. (DEFAULT: 45 Degrees)
%
% 
%roiMask: n x m logical containi
%OUTPUT: 
% subRois: n x m x 3 logical array containing the front, back, and side, subRegional masks 
%          where n and m are the ny and nx of the image and the 3 pages are
%          neurite front, back, and side respectively. These regions are defined 
%          by the local vector ( there is some averaging of this via a user defined  
%          input called vectLength ) at the thickest point  
%          point of a user-selected length from the broad protrusion
%          participating in the longest path- the default here is 20 um and 
%          is set currently in GCAfindNeuriteLength. 
% 
%
ny =  imageSize(1); 
 nx = imageSize(2); 
% %% START 
distTrans = frontPixelInfo(:,2); 
useThickestPoint =1; 
if useThickestPoint == 1;
thickestPt = max(frontPixelInfo(:,2)); 

idxPt =find(distTrans == thickestPt);
else 
    idxPt = length(distTrans);
%[yCenter,xCenter] = ind2sub(imageSize,idxPt); 
end % thickestPoint
% added 20150121 for example for DAC 
makePlot = 1; 
if makePlot == 1
   
    [~,~,distToPlot] = calculateDistance(frontPixelInfo(:,1),[ny,nx],0);
    
    figure; 
    setAxis
    distTransInMic = distTrans.*0.216;
      distTransInMic = distTransInMic(1:end-1)';
      distToPlotInMic = distToPlot.*0.216;
    scatter(distToPlotInMic,distTransInMic,'k','filled'); 
    % mark thickest point 
    scatter(distToPlot(idxPt).*0.216,distTrans(idxPt)*0.216,'r','filled');
   % n = length(distTrans); 
    %distTransInMic = distTrans.*0.216; 
    % for now just find the first radius where the value is equal to 2
    % after the max. 
   
    idxLow =  find(distTransInMic<2 & distToPlot>distToPlot(idxPt),1,'first');     
    line([distToPlot(idxLow).*0.216,distToPlot(idxLow).*0.216],[0,6]); 
   line([0 30],[2,2],'color','k')
    ylabel({'Shortest Distance' ; 'to Veil/ Stem Edge' ; '(um)'});
     xlabel({'Distance Along Neurite Length Path' ;' From Tip of Leading Protrusion ';'(um)'});
    close gcf
 %  sd_spline= csaps(linspace(1,n,n),distTrans,0.01);
%sd=ppval(sd_spline,linspace(1,nTime,nTime));
end     
% get the indices surrounding the vector: length of this is user defined
idx = idxPt-vectLength:idxPt+vectLength;
% make sure to filter out pieces that extend beyond the frontPixelInfo
idx = idx(idx>0); 
idx = idx(idx<=length(frontPixelInfo));

% get slope of local vector 
localVect = frontPixelInfo(idx,1); 
% convert to x,y coordinates 
[yVect,xVect] = ind2sub(imageSize,localVect); 

% make a line along the local vector defined around the thickest point of
% the neurite. 
% slope
mMainAxis = (yVect(1) - yVect(end))./(xVect(1)-xVect(end)); % 1 should hopefully
%be the point on the vector closer to the tip of the protrusion. 
% y-inter
[~,yAll]=meshgrid(1:nx,1:ny);
% calculate the local vector direction back: though note you might end in a
% region where the local direction changes- even still for now just model
% such that goes back. 

deltX = (xVect(1)-xVect(end));
 
if deltX >0 
    dirVectX = 1; 
else 
    dirVectX = -1; 
end
    
   




%% only need if making assymetric subRois
bMainAxis = yVect(end) -  mMainAxis*(xVect(end)); 



%yMainAxis = repmat(mMainAxis.*[1:nx]+bMainAxis,[ny,1]);
% get the value of xy in the OPPOSITE direction of the local vector toward
% the protrusion
% NOTE this might not actually follow the line of the skeleton if the 
% direction in that region changes. Can make accommodations for this at a 
% later time if we need to. 

 % abitrarily shifting slightly away from point by a value of 3 
if abs(mMainAxis) ~= inf; % there can be a small problem when go to infinity in calculating the yVectOpp
xVectOpp = xVect(end)+3*-dirVectX;
    yVectOpp = mMainAxis*(xVectOpp)+bMainAxis;

else % just take it in the opposite direction on the  yAxis  
    xVectOpp = xVect(end); 
    yVectOpp = yVect(end) - (yVect(1)-yVect(end))./abs((yVect(1)-yVect(end)))     ; 
    warning('mMainAxis went to Inf'); 
end 

% convert to pixels 
  pixVectOpp  = sub2ind(imageSize,ceil(yVectOpp),ceil(xVectOpp));
  %pixVect = sub2ind(imageSize,yVect(end),xVect(end)); % get the pixel 
  % at the vector center. 
  
  % this will be used for testing rois the front will contain the forward
  % pointing vector the back will contain the vector in the opposite
  % direction. 
  
  

angleMain = atand(mMainAxis); 


anglePlus = angleMain+angle;
angleMinus = angleMain-angle;

mPlus = tand(anglePlus); 
mMinus = tand(angleMinus); 
% if either slope = infinity just scew it slightly so doesn't error 
if abs(mPlus) == inf
    mPlus = 1e6; % give it a pretty big number for the slope; 
    warning('Slope was equal to infinity : setting slope to 1000000'); 
end 
if abs(mMinus) == inf
    mMinus = 1e6; 
    warning('Slope was equal to infinity : setting slope to 1000000'); 
end 
  bPlus = yVect(end)-mPlus*xVect(end); % maybe want to put the intersection point at the thickest point instead of hte end point...FIX
  bMinus = yVect(end)-mMinus*xVect(end); 
  
yPlus = repmat(mPlus.*(1:nx)+bPlus,[ny,1]); % 
yMinus = repmat(mMinus.*(1:nx)+bMinus,[ny,1]); 
% divide the data into quadrents. 
                dAbove = (yAll<yPlus); %  
                dBelow = (yAll>yMinus); %   
                
           dAboveP =  yAll>yPlus ;    
           dBelowM = yAll<yMinus;     
subRois(:,:,1) = (dAbove & dBelow & roiMask); 
subRois(:,:,2) = (dAboveP & dBelowM &  roiMask);
subRois(:,:,3) = (dAbove & dBelowM & roiMask); 
subRois(:,:,4) = (dAboveP & dBelow & roiMask); 
%subRois(:,:,3) = roiMask - subRois(:,:,1) - subRois(:,:,2);

% for each subRoi assign front middle or side by testing if the point along the  local
% vector lies in the region. the front will be the direction toward the protrusion
% the back will be the opposite direction.. 

 roiBackGC = subRois(:,:,arrayfun(@(i) ~isempty(intersect(find(subRois(:,:,i)==1),pixVectOpp)),1:4)); 
 if isempty(roiBackGC); 
 roiBackGC = zeros([ny,nx]); % just to make sure doesn't error 
 end 

 roiFrontGC = subRois(:,:,arrayfun(@(i) ~isempty(intersect(find(subRois(:,:,i)==1),localVect(1))),1:4)); 
  if isempty(roiFrontGC) 
    roiFrontGC = zeros([ny,nx]); 
 end 
  
 roiSidesGC = roiMask - roiBackGC-roiFrontGC; 

% everything else is the side
%roiSidesGC = roiMask - roiFrontGC - roiBackGC; 

% the only problem with this is it can have extra pieces 
% want only to the first intersect with the body 
% or the CC that includes the central point 
 test = bwconncomp(roiSidesGC); 
 if test.NumObjects > 1; % more than one connected component
% find those CCs that don't intersect with the central point to discard 
toDiscard = vertcat(test.PixelIdxList{cellfun(@(x) isempty(intersect(x,localVect(end))),test.PixelIdxList)}); 
  
 roiSidesGC(toDiscard) = 0; 
 roiBackGC(toDiscard) =1; % typically add it to the back of the subroi
  end 
 % I think for now just store these in an array 
 % page 1 = front 
 % page 2= back 
 % page 3 = sides 
 
subRois(:,:,1) = roiFrontGC; 
subRois(:,:,2) = roiBackGC; 
subRois(:,:,3) = roiSidesGC; 



    

end 



 


%
% 
% 





