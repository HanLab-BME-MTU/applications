function [ subRoisOut,xVect,yVect,pixelInfoGC,defineGCPlot,idxCMap,cmap,distTransInMicGC,params,GCDetect] = GCASubRegionalAutoGrowthConeExtraction(frontPixelInfo,roiMask,varargin)
%GCAgetGrowthConeSubRegions: This function makes a front subregion to the growth cone
% finding the largest distance transform along a line of pixels from the
% veil/stem estimated (calculated in previous step), finding the local
% direction of the skeleton in this region, and adding +/- 45 directional
% degrees to this region in order to automate a subroi.
%
%
%INPUT:

%frontPixels: (REQUIRED) n x 2 double array of pixel n = the number of pixels as defined by
%                         the distCutOff in calculateDistance
%                         column 1 is the pixel indices and column 2 is the distance
%                         transform from the veil/stem estimation (large protrusions)
%
%roiMask:     (REQUIRED) n x m logical containing the roiMask to slice (either the full
%
%
%
%GCFinder:    (PARAM)    structure with fields (DEFAULT empty)
%                                 if not empty will find growth cone and
%                                 give those subRois corresponding to
%                                 growth cone and stem.
%                                 if empty will default to thickest point
%                                 regions for further subRoi partitioning
%
%         GCFindNeckCutOff:     scalar : thickness of the neck of the growth cone
%                                 for the GCFinder cut-off in um (DEFAULT 2 um) - The
%                                 GCFinder will find the first instance of
%                                 this value after the first maximum
%                                 (from the neurite tip) in the distance transform
%                                 along the longest path of the neurite.
%
%         GCFindMinLength:      scalar : the smallest allowed length of a growth cone
%   
%         GCFindThickPt:        scalar: assumes the local max is at least this thick
%
%angle:  (PARAM)             scalar: angle for partitioning
%vectLength:   (PARAM) DEFAULT 4
%
%
%
%OUTPUT:
% subRois: r x c x 3 logical array containing the front, back, and side, subRegional masks
%          where r and c are the ny and nx of the image and the 3 pages are
%          neurite front, back, and side respectively. These regions are defined
%          by the local vector (there is some averaging of this via a user defined
%          input called vectLength) at the thickest point
%          point of a user-selected length from the broad protrusion
%          participating in the longest path- the default here is 20 um and
%          is set currently in GCAfindNeuriteLength.
%
%
% pixelInfoGC:  front pixels only including the GC (for input to further subregion making)
%
% defineGCPlot: a figure handle plotting the distance along the neurite longest
%               path versus the distance transform along this path, with detection
%               cut-off for growth cone marked.
%

ip = inputParser;

ip.CaseSensitive = false;
ip.KeepUnmatched = true;

ip.addRequired('frontPixelInfo',@isnumeric);
ip.addRequired('roiMask',@islogical);

ip.addParameter('GCFinder',true,@islogical);
ip.addParameter('GCFindNeckCutOff',2,@isscalar); % the thicknes cut-off for the veil stem. 
ip.addParameter('GCFindMinLength',5,@isscalar); % the smallest allowed length of a growth cone
ip.addParameter('GCFindThickPt',1.5,@isscalar); % assumes the local max is at least this thick

ip.addParameter('vectLength',4,@isscalar);
ip.addParameter('angle',90,@isscalar);

ip.addParameter('ShowFig','off'); 
ip.addParameter('maxCMapValue',10); % in um


ip.addParameter('maxDistProfile',10); 
 
ip.parse(frontPixelInfo,roiMask,varargin{:});
params = ip.Results; 
%%
GCFinder = ip.Results.GCFinder;
vectLength= ip.Results.vectLength;
angle = ip.Results.angle;

%% Inititiate 
imageSize = size(roiMask);
ny =  imageSize(1);
nx = imageSize(2);

if size(frontPixelInfo,1)==2
distTrans = frontPixelInfo(:,2); % make sure to input as um (or can calculate this internaly as have the roimaks)
end 

defineGCPlot = []; 
idxCMap = []; 
 
%% % useThickestPoint =1;
% if useThickestPoint == 1;
%     thickestPt = max(frontPixelInfo(:,2));
%     
%     idxPt =find(distTrans == thickestPt);
% else
%     idxPt = length(distTrans);
%     %[yCenter,xCenter] = ind2sub(imageSize,idxPt);
% end % thickestPoint
%% Find the first point along the longest path that the distance transform
% hits user specified cut-off value
% FOR THE GC FINDER Fix Later 
% 
%% 

   % Get the distance values along the pixels
[~,~,distToPlot] = calculateDistance(frontPixelInfo(:,1),[ny,nx]); 
% convert to xy coords 
[y,x] = ind2sub([ny,nx],frontPixelInfo(:,1)); 
distTrans = bwdist(~roiMask); 
distTransInMic = distTrans.*0.216; % make a variable
distTransInMic = distTransInMic(frontPixelInfo(:,1)); 
distTransInMic = distTransInMic(1:end-1)';
distToPlotInMic = distToPlot.*0.216; % make a variable

  % initiate the colormap (% make default)
              cmap = jet(128);
              
              % get the average velocity values for all the windows in the current frame
              % and assign a color based on the the mapper.
              plotValues = double(distTransInMic)'; 
              
              mapper=linspace(0,ip.Results.maxCMapValue,128)'; % 
              D=createDistanceMatrix(plotValues,mapper);
              [sD,idxCMap]=sort(abs(D),2);


              if GCFinder;
                  
                  paramNeckCutOff = ip.Results.GCFindNeckCutOff; % in um % MAKE INPUT
                  
                  paramMinGCLength = ip.Results.GCFindMinLength;
                  
                  paramMinGCThickestPt  = ip.Results.GCFindThickPt;
                  
                  % create spline
                  sd_spline= csaps(double(distToPlotInMic),double(distTransInMic));
                  
                  %Find the location of extrema
                  extrema = fnzeros(fnder(sd_spline));
                  
                  extrema = extrema(1,:); % remove the intervals
                  
                  extrema = extrema((extrema ~= ... %These will always be at first or last breaks in spline
                      sd_spline.breaks(1)) & (extrema ~= sd_spline.breaks(end)));
                  sd = ppval(sd_spline,distToPlotInMic);
                  sdExt=ppval(sd_spline,extrema);
                  
                  %  test = fnval(fnder(sd_spline,2),extrema);
                  %Determine whether each extrema is maximum or minimum
                  isMax = fnval(fnder(sd_spline,2),extrema) < 0;
                  % initiate GCLength
                  %GCLength=0;
                  toPlotExt = [extrema(isMax)',sdExt(isMax)'];
                  
                  
                  finalTest = isMax & sdExt>paramMinGCThickestPt ;
                 if sum(finalTest)~=0
                      
                      GCLengthC =0;
                      % iterate until you find a GC length that meets the GC length
                      % requirement
                      while GCLengthC<paramMinGCLength
                          
                          %Find the first maximum that meets GC thickness requirement
                          iBackMax = find(finalTest,1,'first');
                          
                          distFirstMax = extrema(iBackMax) ;
                          
                          % find all lengths of GC that are greater than the first max meeting
                          % the paramGCThickestPt criteria and the paramNeckCutOff
                          % potentialGCLengths = distToPlotInMic(distTransInMic<paramNeckCutOff & distToPlotInMic>distFirstMax);
                          
                          idxNeck = find(distTransInMic<paramNeckCutOff & distToPlotInMic>distFirstMax,1,'first');
                          idxPt = idxNeck; % point for subRegions.
                          % test to make sure meets the paramMinGCLength criteria
                          
                          % make sure the point found is long enough along the distance
                          GCLengthC = distToPlotInMic(idxNeck);
                          finalTest= finalTest(iBackMax+1:end); % take out the previous maximum
                          extrema = extrema(iBackMax+1:end);
                      end
                      GCDetect =true; 
                 else 
                     display('No GC Found');
                     GCDetect =false;
                      
                      % find the first GC length > paramMinGCLength
                      % GCLength = potentialGCLengths(find(potentialGCLengths>paramMinGCLength,1,'first'));
                      %idxNeck = find(distToPlotInMic == GCLength);
                      %idxNeck =  find(distTransInMic<paramNeckCutOff & distToPlotInMic>distFirstMax,1,'first');
                      
                      %     idxPt =find(distTransInMic == thickestPt);
                      idxPt = length(frontPixelInfo)-1; % use the end of the input data
                      pixelInfoGC = frontPixelInfo;
                      idxNeck = idxPt; 
                      
                      
                 end 
                  else % idxPt defaults to the thickest point
                      
                      
                      %     thickestPt = max(distTransInMic);
                      %     idxPt =find(distTransInMic == thickestPt);
                      idxPt = length(frontPixelInfo)-1; % use the end of the input data
                      pixelInfoGC = frontPixelInfo;
                      idxNeck = idxPt;
                
                
              end % if GCfinder
    %idxPt = find(GC == max(GC));
    %idxPt = idxNeck;
%% Make Plot      
    defineGCPlot = setAxis(ip.Results.ShowFig); 
    
    % plot the parmaters for the growth cone defition.
    hold on 
    if ip.Results.GCFinder && GCDetect 
        line([paramMinGCLength,paramMinGCLength], [0 10],'color','k','linewidth',2,'LineStyle','--');
        text(paramMinGCLength+0.1,5,['GC Length Greater Than' num2str(paramMinGCLength)]);
        hold on
        line([0,ip.Results.maxDistProfile], [paramNeckCutOff,paramNeckCutOff],'color','k','linewidth',2,'LineStyle',':');
        line([0,ip.Results.maxDistProfile],[paramMinGCThickestPt,paramMinGCThickestPt],'color','k','linewidth',2,'LineStyle','-');
        text(12.5,paramMinGCThickestPt+0.1,['Local Max : Greater than ' num2str(paramMinGCThickestPt)]);
    end
     
     % Plot Events Selected 
              for iColor = 1:length(cmap)
                 if ~isempty(plotValues(idxCMap(:,1)==iColor));
                 scatter(distToPlotInMic(idxCMap(:,1)==iColor),distTransInMic(idxCMap(:,1)==iColor),50,cmap(iColor,:),'filled'); 
                 end
              end

     %sdExt=ppval(sd_spline,extrema);
     % plot the local maximum found for troubleshooting
     if ip.Results.GCFinder && GCDetect
         scatter(toPlotExt(:,1),toPlotExt(:,2),50,'k','x');
     end
     
     
    scatter(distToPlotInMic(idxNeck),distTransInMic(idxNeck),100,'filled','k');
    if ip.Results.GCFinder
        plot(distToPlotInMic,sd,'k');
    end
    
    
    if     ip.Results.GCFinder
        line([distToPlotInMic(idxNeck),distToPlotInMic(idxNeck)],[0,10],'color','k','Linewidth',4);
    end
    
    scatter(distToPlotInMic(idxNeck),distTransInMic(idxNeck),100,'filled','k','MarkerEdgeColor','w');
    
    if ~GCDetect
      title('No GC Detected');  
    end
    
    ylabel({'Shortest Distance' ; 'to Veil/ Stem Edge' ; '(um)'});
    xlabel({'Distance Along Neurite Length Path' ;' From Tip of Leading Protrusion ';'(um)'});
    axis([0 ip.Results.maxDistProfile 0 10]); 
    angle = 90;
    pixelInfoGC = frontPixelInfo(1:idxNeck,:);
    distTransInMicGC = distTransInMic(:,1:idxNeck)'; 

%% MAKE THE PRIMARY GROWTH CONE MASK

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
%% Fix Problem with Slope if necessary
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
 %dAbove = (yAll<yMinus)
 dBelow = (yAll>yMinus); %

% subRois(:,:,1) = dAbove & roiMask; 
% subRois(:,:,2) = dBelow & roiMask; 

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

roiBack= subRois(:,:,arrayfun(@(i) ~isempty(intersect(find(subRois(:,:,i)==1),pixVectOpp)),1:4));
if isempty(roiBack);
    roiBack= zeros([ny,nx]); % just to make sure doesn't error
end

roiFront = subRois(:,:,arrayfun(@(i) ~isempty(intersect(find(subRois(:,:,i)==1),localVect(1))),1:4));
if isempty(roiFront)
    roiFront = zeros([ny,nx]);
end

%roiSidesGC = roiMask - roiBack-roiFront;

% everything else is the side
roiSidesGC = roiMask - roiBack - roiFront;

% the only problem with this is it can have extra pieces
% want only to the first intersect with the body
% or the CC that includes the central point
test = bwconncomp(roiSidesGC);
if test.NumObjects > 1; % more than one connected component
    % find those CCs that don't intersect with the central point to discard
    toDiscard = vertcat(test.PixelIdxList{cellfun(@(x) isempty(intersect(x,localVect(end))),test.PixelIdxList)});
    
    roiSidesGC(toDiscard) = 0;
    roiBack(toDiscard) =1; % typically add it to the back of the subroi
end
% I think for now just store these in an array
% page 1 = front
% page 2= back
% page 3 = sides

% for now automatically find floaters and just get rid of 
roiFront = getLargestCC(roiFront);

subRoisOut(:,:,1) = roiFront;
subRoisOut(:,:,2) = roiBack;
subRoisOut(:,:,3) = roiSidesGC;
%%
end











