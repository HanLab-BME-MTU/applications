function [ costMat,nonlinkMarker,movieInfo] = costMatFilo( movieInfo,neuriteEdge,neuriteMaskTM1,costMatParams,iFrame,saveDir)

%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%INPUT  movieInfo    : somewhat different from typical movieInfo of single 
%                      particle tracker designed to store more information about the 
%                      detected object: an nFrame structure (movieInfo(iFrame)
%                      with field 
%                           .filoInfo(ifilo) nDetectedFilo structure with
%                           fields:
%                               .pi
%       
%                     
%          
% 
%       costMatParams: Structure with the following fields:
%             .searchRadius: Maximum distance between two features in two
%                          consecutive time points that allows linking 
%      img (should make optional for troubleshooting)     
%OUTPUT costMat      : Cost matrix.
%       noLinkCost   : Cost of linking a feature to nothing, as derived
%                      from the distribution of costs.
%       nonlinkMarker: Value indicating that a link is not allowed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CHECK INPUT


searchRadius = costMatParams.searchRadius;
%% Calculate Linking Costs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Distance matrix calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first just calculate the distance between all mean response fits in the first frame
% and all mean response fits in the second frame: possibly use this for 
% search radius cut-off

%get number of features in the 2 time points
n = numel(movieInfo(1).filoInfo); % num filo detect frame 1 
m = numel(movieInfo(2).filoInfo); % num filo detect frame 2

% extract fit coords for current/Next frame % try using near the base 
% would expect very little movement there as long as the total body NOT 
% moving significantly 

%% THINKING: at first tried to track the end point fits a bit silly 
% as these are the most mobile: let's try to simply track the base of the 
% filo for now: this should be relatively stable unless there is
% significant motion of the neurite forward: if there is this motion is
% something that we could potentially correct for.- so far this does NOT 
% seem to be the case %Test for centroid mobility or movement of the skelton?


% coordsCurr = vertcat(movieInfo(1).filoInfo(:).endpointCoordFitXY);
% coordsNext = vertcat(movieInfo(2).filoInfo(:).endpointCoordFitXY); 
%replicate x,y-coordinates at the 2 time points to get n-by-m matrices
% x1 = repmat(coordsCurr(:,1),1,m);
% y1 = repmat(coordsCurr(:,2),1,m);
% x2 = repmat(coordsNext(:,1)',n,1);
% y2 = repmat(coordsNext(:,2)',n,1);
%% Get the coordinates of the filo at the base of the filopodia 
% for now just use the coordinate at the base of the filo. much more stable
% then can just have costs based on distances as would not expect this to
% fluctuate that much: might fluctuate .  


x1 = repmat(arrayfun(@(x) x.Ext_coordsXY(1,1),movieInfo(1).filoInfo)',1,m);
y1 = repmat(arrayfun(@(x) x.Ext_coordsXY(1,2),movieInfo(1).filoInfo)',1,m);
x2 = repmat(arrayfun(@(x) x.Ext_coordsXY(1,1),movieInfo(2).filoInfo),n,1); 
y2 = repmat(arrayfun(@(x) x.Ext_coordsXY(1,2),movieInfo(2).filoInfo),n,1); 

%%
if costMatParams.predict == 1; 


    
% project the xy coordinates forward based on the boundary evolution 
neurite1 = neuriteEdge{1}; % in the end make part of movie info. 
neurite2 = neuriteEdge{2}; 


% try contour
%  currOutline = contourc(double(currMask),[0 0]);
%     currOutline = separateContours(currOutline);%Post-processing of contourc output
%     currOutline = cleanUpContours(currOutline);    
%     currOutline = currOutline{1}';%We know we only have one object...


% maybe just load img here for now 
% img = double(imread(movieInfo(1).imgPointer)); 
% neurite1Mask = zeros(size(img)); 
% idx = sub2ind(size(img),neurite1(:,1),neurite1(:,2)); 
% neurite1Mask(idx) = 1; 
% neurite1Mask = imfill(neurite1Mask); 
%idx2 = sub2ind(size(img),yAtBase,xAtBase);


% 
%maskBeforeConn([idx;idx2']) = 1; 
    



inputParam.TOLERANCE = 30; 
inputParam.ISCLOSE = 0; 
inputParam.batch_processing = 1; 
inputParam.dl_rate = 30; 
inputParam.time_idx = 2; 
inputParam.pixel_list_last = [neurite1(:,2) neurite1(:,1)] ; 
inputParam.pixel_list = [neurite2(:,2),neurite2(:,1)] ; % xy format 
inputParam.maskTM1 = neuriteMaskTM1; 
%inputParam.maskTM1 = imfill(neurite1)     ; % this was previously a bug...figure out why it wasn't a problem. 04/07 
[output] = prSamProtrusion(inputParam); 

% take the end point coords and predict the new coords 
xAtBase = arrayfun(@(x) x.Ext_coordsXY(1,1),movieInfo(1).filoInfo); 
yAtBase = arrayfun(@(x) x.Ext_coordsXY(1,2),movieInfo(1).filoInfo); 



xf = neurite1(:,2);
yf = neurite1(:,1);


stre = [1 1 1; 1 1 1; 1 1 1]; 
    % for now just loop through as need to keep index... 
for iBase = 1:length(xAtBase)
    
    test = zeros(size(img)); 
    idxBase = sub2ind(size(img),yAtBase(iBase),xAtBase(iBase)); 
    test(idxBase) = 1; 
    test = imdilate(test,stre); 
    % find those that overlap with edge
    idxDil = find(test==1); 
    idxEdge =  sub2ind(size(img),yf,xf);
    [idxToPropForTemp,~,idxListTemp]= intersect(idxDil,idxEdge); 
    if ~isempty(idxToPropForTemp)  % currently can be empty at this point 
        % because I don't always deal adequately with branch points... 
        
    idxToPropForFinal(iBase) = idxToPropForTemp(1); % if mor than one just take the first
    idxListFinal(iBase) = idxListTemp(1); 
    else 
        idxToPropForFinal(iBase) = NaN; 
        idxListFinal(iBase) = NaN; 
    end 
end 



% filo coords on neurite edge to propagate forward:
[yProp, xProp] = ind2sub(size(img),idxToPropForFinal);



%
idxNoPred = find(isnan(yProp));
idxPred = find(~isnan(yProp));

% get just the translations for the points want to project forward
%   transY = output.translate_output(idxListFinal(~isnan(idxListFinal)),2); % switched? xy 
% transX = output.translate_output(idxListFinal(~isnan(idxListFinal)),1);

 transY = output.translate_output(:,2); 
  transX = output.translate_output(:,1); 


xPredAll = neurite1(:,2) +transX; 
yPredAll = neurite1(:,1) +transY; 

% xPredAll = xProp(~isnan(xProp))' + transX;
% yPredAll = yProp(~isnan(yProp))' + transY;
% 

% get back the nans for indexing
% xFinal = nan(length(m))';
% yFinal = nan(length(m))';
% 
% % if
% for iPred = 1:length(xPredAll)
%     xFinal(idxPred(iPred)) = xPredAll(iPred);
%     yFinal(idxPred(iPred)) = yPredAll(iPred);
% end
% 
% % if no prediction because end point too far away from neurite edge
% % just use 
% for iNoPred = 1:length(idxNoPred)
%     xFinal(idxNoPred(iNoPred)) = xAtBase(idxNoPred(iNoPred));
%     yFinal(idxNoPred(iNoPred)) = yAtBase(idxNoPred(iNoPred));
% end 

% put all prediction and necessary non-predictions in a matrix to do the 
% dist calc
% xPredFinal = repmat(xFinal',1,m); 
% yPredFinal = repmat(yFinal',1,m);
% x1 = xPredFinal; % change the location of that in the first frame to the prediction
% y1 = yPredFinal; % 


% save the information in the filo Info for now so can extract 
% for iFilo = 1:numel(filoInfo)
%     movieInfo(1).filoInfo(iFilo).predictedDisp = [transX,transY];  
% end

end % if costMatParams.Predict

%%
%% sanity check!! 
troubleshoot =0;

if troubleshoot == 1; 
    
    % first plot all 
    [ny,nx] = size(img); 
    setFigure(nx,ny); 
    imshow(img,[]) ; 
    hold on 
    
    for i = 1:length(neurite1(:,1))
     line([neurite1(i,2) xPredAll(i)],[neurite1(i,1) yPredAll(i)],'Color','g')
    end
 
    
    
    
    
%    [ny,nx] = size(img); 
% setFigure(  nx , ny); 
% imshow(img,[]) 
 hold on 
 plot(neuriteEdge{1}(:,2),neuriteEdge{1}(:,1),'b');
 plot(neuriteEdge{2}(:,2),neuriteEdge{2}(:,1),'y','Linestyle','--'); 
% %scatter(xPredAll,yPredAll,'g'); 
% scatter(xProp,yProp,'c','filled'); 
% 
% xAtBase2 = arrayfun(@(x) x.Ext_coordsXY(1,1),movieInfo(2).filoInfo); 
% yAtBase2 = arrayfun(@(x) x.Ext_coordsXY(1,2),movieInfo(2).filoInfo); 
% % 
% scatter(xAtBase2,yAtBase2,'y','filled'); 
% scatter(xPredAll,yPredAll,'g'); 
% if no predication the distance needs to be from the original
% coords...fix.. should not be a problem if I get the damn branching taken 
% care of. if x : nana
% xPred = repmat(xPredAll,1,m);
%yPred= repmat(yPredAll,1,m);



% xPropOnEdge = xProp(~isnan(xProp))'; 
% yPropOnEdge = yProp(~isnan(yProp))'; 
% 
% 
% 
% for i = 1:length(xPropOnEdge)
%     line([xPropOnEdge(i) xPredAll(i)],[yPropOnEdge(i) yPredAll(i)],'Color','g')
% end 
 if ~isdir([saveDir filesep 'troubleshoot']) 
     mkdir([saveDir filesep 'troubleshoot']) 
 end 
%saveas(gcf,[saveDir filesep 'troubleshoot' filesep  num2str(iFrame) '.esp'],'psc2'); 
saveas(gcf,[saveDir filesep 'troubleshoot' filesep num2str(iFrame) '.png']); 
save([saveDir filesep 'troubleshoot' filesep 'output_Frame' num2str(iFrame) '.mat'],'output'); 

close gcf
end
%%

% find potential candidates based on search radius from predicted base
% position
%calculate the square distances between features in time points t and t+1
distMat = (x1-x2).^2 + (y1-y2).^2;

%assign NaN to all pairs that are separated by a distance > searchRadius
distMat(distMat>searchRadius^2) = NaN;


% for now just normalize here to make on order with the angle term 
distMat = distMat./max(distMat(:)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate angle between the two filo vectors at the base of the filopodia
% for cost param: larger angle between two vectors = higher cost 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

costParams = nan(n,m,3); % initiate # filo frame 1 x # filo frame 2 array for cost matrices

% get the IDs of all potential links 
[IDCurrent,IDNext] = find(~isnan(distMat)); 

% for each particle find potential candidates in next frame and calculate
% cost (not sure this is quickest way to do this in loops, likely better to 
% do all at once but as I have to to loop to make this calc anyway just do here
% )... maybe will do if decide
% maybe decide to change later for elegence. 

for ifilo = 1:n % for all filo in current frame
    
    % get all candidates that fall within a specific range
    IDNextOK = IDNext(IDCurrent==ifilo);
    
    % get relavant info and calculate new cost matrice
    %x1 = movieInfo(1).xCoord(ifilo,1);
    %y1 = movieInfo(1).yCoord(ifilo,2);
    
    fit = round(movieInfo(1).filoInfo(ifilo).Ext_length); % get the fit Length (in pixels)
   % if ~isempty(fit)
   if fit ~=0; % added 20141129 think changed the output of a zero length from [] to 0 
    % divide the fit length by two to get the length value of the filo 'base'
    % to be used for angle calculation. calc the delta y/delta x and mag
    % describing the filo vector
    
    deltaXBaseCurr=  movieInfo(1).filoInfo(ifilo).Ext_coordsXY(round(fit/2),1) - movieInfo(1).filoInfo(ifilo).Ext_coordsXY(1,1);
    deltaYBaseCurr = movieInfo(1).filoInfo(ifilo).Ext_coordsXY(round(fit/2),2)- movieInfo(1).filoInfo(ifilo).Ext_coordsXY(1,2) ;
    magBaseCurr = sqrt(deltaXBaseCurr^2 +  deltaYBaseCurr^2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % troubleshoot plots
%         % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         troubleshoot = 0;
%         if troubleshoot == 1 
%            
%             
%             currentFilo = movieInfo(1).filoInfo(ifilo).Ext_coordsXY; 
%             
%             plot(currentFilo(:,1),currentFilo(:,2),'g');   
%             
%             hold on
%         end 
    
    % again non-elegant way to do it but just run with it
    % find each candidate link for given filo and calculate the various costs 
    % associated with that linkage
   
        
    for iCan = 1:length(IDNextOK)
        
        % get the idx of the candidatate from list
        idxCan = IDNextOK(iCan);
         fit = round(movieInfo(2).filoInfo(idxCan).Ext_length);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate Cos Angle between current frame and viable candidates
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % divide the fit length by two to get the length value of the filo 'base'
        % to be used for angle calculation. calc the delta y/delta x and mag
        % describing the filo vector in the 2nd frame
      %  if ~isempty(fit) 
         if fit~=0 % changed 20141129
           
            deltaXBaseNext =  movieInfo(2).filoInfo(idxCan).Ext_coordsXY(round(fit/2),1) - movieInfo(2).filoInfo(idxCan).Ext_coordsXY(1,1);
            deltaYBaseNext = movieInfo(2).filoInfo(idxCan).Ext_coordsXY(round(fit/2),2)- movieInfo(2).filoInfo(idxCan).Ext_coordsXY(1,2) ;
            magBaseNext = sqrt(deltaXBaseNext.^2 +  deltaYBaseNext.^2);
            
            % calculate angle between the two filo vectors
            cos_ATBase = (deltaXBaseCurr*deltaXBaseNext + deltaYBaseCurr*deltaYBaseNext)/(magBaseNext*magBaseCurr);
            
            % save as a costParam
            costParams(ifilo,idxCan,1) = cos_ATBase;
            
            if troubleshoot == 1
                %plot the candidate filo
                canFilo = movieInfo(2).filoInfo(idxCan).Ext_coordsXY;
                plot(canFilo(:,1),canFilo(:,2),'y');
                
                % plot the search radius
                [x,y] = cylinder(searchRadius);
                coord1 = movieInfo(1).filoInfo(ifilo).Ext_coordsXY(1,:);
                plot(x+coord1(1,1),y+coord1(1,2),'r');
                
                
            end % if troubleshoot
        
        else 
            % eventually filter out the ones that have such a poor fit. 
            % however for now just filter here. 
            costParams(ifilo,idxCan,1)= NaN ; 
        end  
       
    end % for iCan
    else % go on to next
    end % isempty
end % for iFilo
%% 
costMat = distMat + (1-costParams(:,:,1)); 
%% Birth and death

%append matrix to allow birth and death
% jonas, 10/09: fix for non-sparse tracker
% if isstruct(prevCost)
%     prevCostMax = prevCost.max;
% else
%     prevCostMax = max(prevCost(:));
% end
% if ~isnan(prevCostMax) && prevCostMax ~= 0
%     maxCost = 1.05*prevCostMax;
% else
% max cost is determined by taking the 80th percentile of all costs
    maxCost = max(prctile(costMat(:),80));
% end
deathCost = maxCost * ones(n,1);
birthCost = maxCost * ones(m,1);

%generate upper right and lower left block
deathBlock = diag(deathCost); %upper right
deathBlock(deathBlock==0) = NaN;
birthBlock = diag(birthCost); %lower left
birthBlock(birthBlock==0) = NaN;

%get the cost for the lower right block
costLR = min(min(min(costMat))-1,-1);
lrBlock = costMat';
lrBlock(~isnan(lrBlock)) = costLR;

%append cost matrix
costMat = [costMat deathBlock; birthBlock lrBlock];

%% nonLinkMarker

%determine the nonlinkMarker
nonlinkMarker = min(floor(min(min(costMat)))-5,-5);

%replace NaN, indicating pairs that cannot be linked, with nonlinkMarker
costMat(isnan(costMat)) = nonlinkMarker;

end


 
   
        
        
        




