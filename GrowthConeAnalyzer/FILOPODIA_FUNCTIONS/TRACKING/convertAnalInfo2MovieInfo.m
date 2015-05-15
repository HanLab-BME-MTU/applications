function [movieInfo , analInfoFilt] = convertAnalInfo2MovieInfo( analInfo )
%just convert analInfo to movieInfo for now for historical reasons might want to change
%structures later to make pretty 

analInfoFilt = analInfo; 
for iFrame = 1:numel(analInfo) -1 % for most of these I don't record the last frame make a task to remind myself why

 
filoInfo = analInfo(iFrame).filoInfo; 


     
% here I fix the bad fits to just use the original length 
% NOTE should just save it like this in the fitting function!!! FIX ! 
% fix lengths
 idxBadFit =  find(isnan(vertcat(filoInfo(:).Ext_length)));
    % just use original length 
    lengthsOrig = arrayfun(@(x) length(x.Ext_pixIndicesBack),filoInfo); 
    lengthsOrig = lengthsOrig'; 
    coord = arrayfun(@(x) [x.Ext_endpointCoord(:,2) x.Ext_endpointCoord(:,1)],filoInfo,'uniformoutput',0); 
    for i = 1:length(idxBadFit)
    filoInfo(idxBadFit(i)).Ext_length = lengthsOrig(idxBadFit(i)); 
    filoInfo(idxBadFit(i)).Ext_endpointCoordFitXY = coord{idxBadFit(i)}; 
    end 
    noBranchSteps =1; 
    if noBranchSteps ==1 ; 
        
        
     filoInfoFilt = filoInfo; % just track everything separately %%% NOTE
     %NEED TO MAKE MORE EXPLICIT INPUT 
%        types = vertcat(filoInfo(:).type); 
%        filoInfoFilt = filoInfo(types==0|types==1); 
%      
     
    else 
    
    % filter based on type (just take neurite body based ones to start)
  types = vertcat(filoInfo(:).type); 
  % as a quick fix... make so that when vertcat you don't lose these
  % indices 
  idxFix = find(types==0)  ;
  for i = 1:length(idxFix)
      filoInfo(idxFix(i)).conXYCoords = [NaN NaN NaN]; % NOTE should have this when start the 
  end 
  
  
  filoInfoFilt = vertcat(filoInfo(types <=2));  % keep only body bodystems and lower order branches that might 
  % be close to the neurite body 
  conXYDist = vertcat(filoInfoFilt(:).conXYCoords);
  types2 = vertcat(filoInfoFilt(:).type); 
  % don't include far branches yet 
  filoInfoFilt= filoInfoFilt(~(conXYDist(:,3)>3 & types2==2)); % keep only those secondary branches 
  % that are have a branchpoint 3 pixels 
    end 
  %types2 = vertcat(filoInfoFilt(:).type); 
    % get all main branches and calculate the connectivity index 
    troubleshoot=0; 
   if troubleshoot ==1 
      img = double(imread(analInfo(iFrame).imgPointer)); 
      figure; 
      imshow(img,[]) 
      hold on 
      test = vertcat(filoInfoFilt(:).Ext_pixIndicesBack); 
      mask = zeros(size(img)); 
      mask(test) =1; 
      spy(mask); 
   end 
lengths = vertcat(filoInfoFilt(:).Ext_length);
std = vertcat(filoInfoFilt(:).Ext_std); 
xyCoords= vertcat(filoInfoFilt(:).Ext_endpointCoordFitXY); % be careful do I even need this? I am tracking the base points... 
% % could put the base coords in here instead... but want some orientation
% % information to be built into the cost function so not necessarily
% % straight forward could maybe have that information here... 

% added 20140928
idxNaN = find(arrayfun(@(x) isempty(filoInfoFilt(x).windowIdx),1:length(filoInfoFilt)));
if ~isempty(idxNaN)
% for now fix empty with NaN
for i=1:length(idxNaN); 
    filoInfoFilt(idxNaN(i)).windowIdx = NaN; 

end 
end 
windowNum =arrayfun(@(x) filoInfoFilt(x).windowIdx(1),1:length(filoInfoFilt)); 
%localVeil = arrayfun(@(x) filoInfoFilt(x).localVeil,1:length(filoInfoFilt));  
add = 0.5*ones(length(xyCoords),1); 
movieInfo(iFrame).xCoord = [xyCoords(:,1) add]; 
movieInfo(iFrame).yCoord = [xyCoords(:,2) add]; 
%movieInfo(iFrame).amp = [lengths std]; % remind me why I need movieInfo format again... should try to eradicate 
movieInfo(iFrame).amp = [windowNum' std];
%movieInfo(iFrame).amp = [localVeil' std]; 
 analInfoFilt(iFrame).filoInfo = filoInfoFilt; 



end



end

