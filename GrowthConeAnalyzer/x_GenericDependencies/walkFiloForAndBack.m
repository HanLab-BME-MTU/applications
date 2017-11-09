function [filoInfo] = walkFiloForAndBack(filoInfo,verticesEP,edgePathCoord,maxTh,maxRes,img,int,n_pixFor,bodyMask)
%short function that extracts and the records the response and intensity 
% values of the edgePaths and propogates the coordinates of the filo 
% forward... prepares filoInfo so that it can be input into the linescan
% function to fit endpoint coordinate.
% filoInfo: you can begin to add fields to the structure. might want to
% initiate in the beginning? 
% int: flag to tell what field names to use in the structure.  1 is
% internal else the filo are assumed to be external filopodia 
% n_pixFor: the number of pixels you want to project the filo forward for
% the fit. 2013_07_14 it seems like this is helpful to use different values
% for internal and external... internal might want to be a bit more
% permissive. 
% adds to a filoInfo structure
% if isempty(filoInfo)
%     numFiloPrev=0; 
% else 
% numFiloPrev = numel(filoInfo);
% end
 [ny nx] = size(maxTh) ;


    %% GET AVG RESPONSES AND INTENSITIES FOR FILO COORDS AND WALK FORWARD BASED ON STEERABLE FILTER DIRECTION

% filter for only filopodia of a certain length when to do this?% 
% for now just indicate short filo. 
filterSize = 0; 
shortFilo = zeros(numel(edgePathCoord),1); 
[numFilo] = size(verticesEP); 
for ifilo = 1:numFilo
 %  filoCount = numFiloPrev +ifilo; 
  filoCount = ifilo; 
    filoCoords= edgePathCoord{ifilo};
    if length(filoCoords) > filterSize == 1 
        
    pixIdxBack = sub2ind(size(maxTh),filoCoords(:,1),filoCoords(:,2)); 
    
    pixIdxBack=  flipdim(pixIdxBack,1); 
    filoCoords = flipdim(filoCoords,1);
%%%%%%%%%%%%%% AVERAGE RESPONSE BACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    avgResBack = nan(length(pixIdxBack),1); 
    avgIntBack = nan(length(pixIdxBack),1); 
  
  [yMax,xMax] = size(maxTh);
  
  
    
    
    indicesForMaskBack = zeros(length(pixIdxBack),5); 
  for i= 1:length(pixIdxBack) 
   coordsPer = zeros(5,2); 
   dx = round(cos(maxTh(pixIdxBack(i))));  
   dy = round(sin(maxTh(pixIdxBack(i)))); 
   
   coordsPer(1,1) = filoCoords(i,1); %
   coordsPer(1,2) = filoCoords(i,2); 
   
   
   coordsPer(2,1) = coordsPer(1,1) + dy; 
   coordsPer(2,2) = coordsPer(1,2) + dx; 
   
   coordsPer(3,1) = coordsPer(1,1) -dy; 
   coordsPer(3,2) = coordsPer(1,2) - dx; 
   
   coordsPer(4,1) = coordsPer(2,1) +dy; 
   coordsPer(4,2) = coordsPer(2,2) +dx; 
   coordsPer(5,1) = coordsPer(3,1) - dy; 
   coordsPer(5,2) = coordsPer(3,2) - dx; 
      
     
    idxOutOfBounds = find(coordsPer(:,1) > yMax | coordsPer(:,2) > xMax |...
        coordsPer(:,1) <= 0 | coordsPer(:,2) <= 0);  
    if ~isempty(idxOutOfBounds)
    coordsPer(idxOutOfBounds,:) = nan; 
    end 
    
   
    indicesForAvgBack = sub2ind(size(img),coordsPer(:,1),coordsPer(:,2));
     indicesForMaskBack(i,:) = indicesForAvgBack'; 
    indicesForAvgBack = indicesForAvgBack(~isnan(indicesForAvgBack)); 
    
   
    
    avgResBack(i) = nanmean(maxRes(indicesForAvgBack)); 
    avgIntBack(i) = nanmean(img(indicesForAvgBack));
    
    
 end      
       
   
    
    
    
    
    % NOTE 
    % currently just average the response generically over local environ
    % ie all 8 nn pixels might want to modify in the future.  
    
    % hmmm.. there has to be a function that does this right? 
    % gets all the local area pixels? lets do it the hard way
    
%     for i = 1:length(pixIdxBack) 
%          dist = zeros(size(img)); 
%          dist(pixIdxBack(i)) =1; 
%          dist = bwdist(dist,'chessboard');
%          nn = find(dist == 1); 
%          
%          avgResBack(i) = mean(maxRes([pixIdxBack(i);nn])); 
%          avgIntBack(i) = mean(img([pixIdxBack(i);nn])); 
%         
%     end 
    


         

    %%% WALK FORWARD %%%
    
    
    % project the filopodia forward based on the direction of the response
    % from the steerable filter
    
   n =  n_pixFor ; % how many pixels you will project forward
    
    % initialize (from endpoint) 
    coordsFor = nan(n+1,2);
    dirFeat = nan(n+1,1); 
    
    coordsFor(1,1) = verticesEP(ifilo,1); % y
    coordsFor(1,2) = verticesEP(ifilo,2); % x
    
    % NOTE direction stored in  is perpendicular to filo length
     dirFeat(1) = maxTh(coordsFor(1,1),coordsFor(1,2));
    % find which direction oriented 
    
   
    
      % Define the direction forward from tip 
        filoDeltX = coordsFor(1,2)-filoCoords(end,2);
        filoDeltY = coordsFor(1,1)-filoCoords(end,1); 
        filoMag = sqrt(filoDeltX^2 + filoDeltY^2); 
        
        dX1 = round(cos(dirFeat(1)-pi/2));
        dY1 = round(sin(dirFeat(1)-pi/2)); 
        
        dX2 = round(cos(dirFeat(1)+pi/2)); 
        dY2 = round(sin(dirFeat(1)+pi/2)); 
        
        cosV1Filo = (filoDeltX.*dX1 + filoDeltY.*dY1)/filoMag; 
        cosV2Filo = (filoDeltX.*dX2 + filoDeltY.*dY2)/filoMag; 
        
       test1 = (1-cosV1Filo);
       %test2 =(1-cosV2Filo); 
       
       if abs(test1) < 1
           dX = dX1; 
           dY = dY1; 
           add = -pi/2;
       else 
           dX = dX2; 
           dY = dY2; 
           add = pi/2; 
       end 
       
      
        
       % filoInfo(ifilo).endPointCoords = [verticesEP(ifilo,1),verticesEP(ifilo,2)]; 
        %test = atan(verticesEP(ifilo,1)-yBack(1),verticesEP(ifilo,2)-xBack(1));
        %if test > 0
         %   add = -pi/2 ; else add = +pi/2; 
        %end 
        
       
       coordsFor(2,2) = coordsFor(1,2) + dX; 
       coordsFor(2,1) = coordsFor(1,1) + dY;     
        
        propBasedOnResp = 0; % option to propogate forward based on the 
        % of the response at each pixel (might either be better or get
        % noisy) 

        % use the direction at tip for all projected values imshow(img,[]);

    for i = 2:n 
         if (round(coordsFor(i,1))<=0 || round(coordsFor(i,2))<=0 || round(coordsFor(i,1))>ny || round(coordsFor(i,2))>nx)% 
             % break if hit end of image 
             
            idx = ~isnan(coordsFor(:,1)); % shorten the coords 
            coordsFor = coordsFor(idx,:);
            coordsFor = coordsFor(1:end-1,:); % delete the problem coord
             
             
            break
           
         else 
         if propBasedOnResp == 1
             % recalculate dir based on new coord 
             
                 dirFeat(i) = maxTh(round(coordsFor(i-1,1)),round(coordsFor(i-1,2)));  
                 dX = round(cos(dirFeat(i)+add)); 
                 dY = round(sin(dirFeat(i)+add)); 
                 coordsFor(i+1,2) = coordsFor(i,2) + dX; 
                 coordsFor(i+1,1) = coordsFor(i,1) + dY; 
         else  % project coordinates only based on the direction of the 
             % response at the tip of the filo
                coordsFor(i+1,2) = coordsFor(i,2) + dX; 
                coordsFor(i+1,1) = coordsFor(i,1) + dY; 
         end 
         end 
    end 
     
    idxOutOfBounds = find(coordsFor(:,1) > yMax | coordsFor(:,2) > xMax |...
        coordsFor(:,1) <= 0 | coordsFor(:,2) <= 0);  
    if ~isempty(idxOutOfBounds)
    coordsFor(idxOutOfBounds,:) = nan ; 
    end 
     
  
   pixIdxFor =  sub2ind(size(maxTh),round(coordsFor(:,1)),round(coordsFor(:,2))); 
  
%    %%%% AVERAGE RESPONSE AND INTENSITY FORWARD %%%
%    avgResFor = nan(n+1,1); 
%    avgIntFor = nan(n+1,1); 
%    
%     for i = 1:length(pixIdxFor) 
%          dist = zeros(size(img)); 
%          dist(pixIdxFor(i)) =1; 
%          dist = bwdist(dist,'chessboard');
%          nn = find(dist == 1); 
%          
%          avgResFor(i,1) = mean(maxRes([pixIdxFor(i);nn])); 
%          avgIntFor(i,1) = mean(img([pixIdxFor(i);nn])); 
%     end 



% get average response 2 pixels on each sideof coordinate
% for now just do direction feature of tip %% 
% SHOULD CHANGE TO INCORPORATE INFORMATION OF SCALES!!!
avgResFor = nan(n+1,1); 
avgIntFor = nan(n+1,1); 
% stupid way to do this think about changing. 
coordsPer = zeros(5,2);

indicesForMaskFor = zeros(length(pixIdxFor),5); 

for i= 1:length(pixIdxFor) 

   dx = round(cos(dirFeat(1))); 
   dy = round(sin(dirFeat(1))); 
   
   coordsPer(1,1) = coordsFor(i,1);% maybe add a third d here: therefore save as 
   coordsPer(1,2) = coordsFor(i,2); % first d is the center coord (xy) - maybe just make  coordsPer(i,1,1) etc 
   % where i is the ID along the filo - though the easiest way to do it for
   % visualization is just dilation....also could just get the dilataion
   % along each point and weight center high surrounding lower - in this
   % manner front and behind will be taken in the signal ... could likely
  
   
   
   coordsPer(2,1) = coordsPer(1,1) + dy; % first set of flanking coords 
   coordsPer(2,2) = coordsPer(1,2) + dx; 
   
   coordsPer(3,1) = coordsPer(1,1) -dy; 
   coordsPer(3,2) = coordsPer(1,2) - dx; 
   
   coordsPer(4,1) = coordsPer(2,1) +dy; % second set of flanking coords
   coordsPer(4,2) = coordsPer(2,2) +dx; 
   coordsPer(5,1) = coordsPer(3,1) - dy; 
   coordsPer(5,2) = coordsPer(3,2) - dx; 
     
    idxOutOfBounds = find(coordsPer(:,1) > yMax | coordsPer(:,2) > xMax |...
        coordsPer(:,1) <= 0 | coordsPer(:,2) <= 0);  
    if ~isempty(idxOutOfBounds)
    coordsPer(idxOutOfBounds,:) = nan ; 
    end 
    
    
    
    indicesForAvg = sub2ind(size(img),coordsPer(:,1),coordsPer(:,2));  
    indicesForMaskFor(i,:)  = indicesForAvg';
    indicesForAvg = indicesForAvg(~isnan(indicesForAvg)); 
    
    
    avgResFor(i) = nanmean(maxRes(indicesForAvg)); 
    avgIntFor(i) = nanmean(img(indicesForAvg)); 
 end      
     
  

  %%% SAVE RESP AND INTENSITY INFO %%%
    
    avgResFilo = [avgResBack;avgResFor] ;
    avgIntFilo = [avgIntBack;avgIntFor]; 
   %% Make sure to filter out pixels that extend beyond the body Mask (however need to make sure pixIdxFor is in same order 
   if int ==1 % you need to filter to make sure the pixels don't excede the body estimation
       pixelsForBody = find(bodyMask==1);
       pixIdxFor = intersect(pixIdxFor,pixelsForBody,'stable');
       
       % Addition Added 20151020 : small check to make sure your forward
       % extension is not in the realm of the boundary pixels so you do not 
       % get a false fit. 
       
       yxBound = bwboundaries(bodyMask);
       
       boundaryPixels = sub2ind([ny,nx],yxBound{1}(:,1),yxBound{1}(:,2));
       
       boundMask = zeros([ny,nx]);
       boundMask(boundaryPixels) = 1;
       boundDilate = imdilate(boundMask,strel('disk',4));
       % create a logical index array such that any points within the
       % final fit region are deleted.
       pixDilate = find(boundDilate);
       
       toDelete = [length(pixIdxFor)-3: length(pixIdxFor)];
       [~,idx] = intersect(pixIdxFor,pixDilate); % this will be all
       
       toDelete = intersect(idx,toDelete);
       if ~isempty(toDelete)
           pixIdxFor(toDelete)=[];
       end
   end
    pixIndices = [pixIdxBack; pixIdxFor];
    else % keepinfo about these nubs in the structure for now..may be useful later if they are sites of nascent filopodia
        avgResFilo = nan;
        avgIntFilo = nan;
        pixIndices = nan; 
        dirFeat = maxTh(verticesEP(1,1),verticesEP(1,2)); 
        pixIdxFor = nan; 
        pixIdxBack = nan; 
        shortFilo(ifilo,1) = 1; % save ID to filter our filoInfo later
        indicesForMaskBack = nan(1,5);
        indicesForMaskFor = nan(1,5); 
    end  % isempty filoCoords      
   
    if int == 1;
    toAdd = 'Int_';
    else toAdd = 'Ext_'; 
    end 
    
    % think I want to save everything then filter to fit. 
    filoInfo(filoCount).([toAdd 'response' ]) = avgResFilo; 
    filoInfo(filoCount).([toAdd 'intensities' ]) = avgIntFilo; 
    filoInfo(filoCount).([toAdd 'maskIndices']) = [indicesForMaskBack; indicesForMaskFor];
    %filoInfo(filoCount).([toAdd 'maskPostFit']) =   ; 
    filoInfo(filoCount).([toAdd 'dirAtTip' ]) = dirFeat(:);
    filoInfo(filoCount).([toAdd 'pixIndices' ]) =  pixIndices; % both forward and back
     [y,x] = ind2sub(size(img),pixIndices);
    filoInfo(filoCount).([toAdd 'coordsXY' ]) =[x,y]; % all potential filoCoords;
     
    filoInfo(filoCount).([toAdd 'endpointCoord' ]) = verticesEP(ifilo,:); % based on response
    filoInfo(filoCount).([toAdd 'pixIndicesFor' ]) = pixIdxFor; % for plotting
    filoInfo(filoCount).([toAdd 'pixIndicesBack' ]) = pixIdxBack; % for plotting
    filoInfo(filoCount).([toAdd 'vectFilo']) = [x(end,1)-x(1,1),y(end,1)-y(1,1)]; 
%     filoInfo(filoCount).type = type; 
%     filoInfo(filoCount).numTrackObj = filoCount; % for now these are singles so part of same tracking ob
%   
    
     
    
    
  
end %ifilo
%% filter short filopodia % maybe take out 
% filoInfo = filoInfo(~logical(shortFilo));  % analInfoFilt(iFrame)
% %ID = 1:length(verticesEP) ; 
% %longFilo = ID(shortFilo==0) ; 
% % for now filter by short filo
% % sort by short filo 
%    % analInfo(iFrame).filoInfo = analInfo(iFrame).filoInfo(~logical(shortFilo')); 
% % filter filoInfo to only those that meet the "long filo" criteria 
% 
% pixIdxAllProj = vertcat(filoInfo.pixIndicesFor); 
% pixIdxAllProj= pixIdxAllProj(~isnan(pixIdxAllProj)); 
% testProj = zeros(size(img)); 
% testProj(pixIdxAllProj) = 1; 
% 
% 
% pixIdxAllBack = vertcat(filoInfo.pixIndicesBack); 
% pixIdxAllBack = pixIdxAllBack(~isnan(pixIdxAllBack)); 
% testBack = zeros(size(img)); 
% testBack(pixIdxAllBack) = 1; 
% 
% filterVert = verticesEP(~logical(shortFilo),:); 


end   


