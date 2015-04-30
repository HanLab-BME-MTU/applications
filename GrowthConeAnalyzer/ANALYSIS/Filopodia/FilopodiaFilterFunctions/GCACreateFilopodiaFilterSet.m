function [ filoFilterSet,filterParams] = GCACreateFilopodiaFilterSet(analInfo,filterType)
%GCACreateFilopodiaFilterSet:
% create a logical filterset for filopodia analysis - ... maybe just have the
% input be a filoFilterSet created previously (with filter types and
% associated analysis) 
%
%
% INPUT:
%      analInfo:     REQUIRED:
%
%
%      filterType:   PARAM:    char   Currently Either 'ConnectVeil_LengthInt', 'ConnectVeil_DensityOrient','Branch',
%                                     'Internal' , 'ManualDefine', 
%
%                                     'VeilConnect': Will implement the Set
%                                      Defaults for extracting the veil
%                                      parameters
%
% OUTPUT:
%       filoFilterSet:   rx1 cell of logical filters for the filopodia for each frame
%                        where r1 is the number of frames 
%                        
%                             
%
%% SET THE PARAMS BASED ON INPUT CHOICE
switch filterType
    case 'ConnectToVeil_LengthInt';
        filterParams.filoTypes = [0,1]; % 0 order attached to veil (no Branch), 1st order attached to a veil with a branch
        filterParams.filterByFit = 1; % turn on filtering by filopodia fit.
        filterParams.filterByBundleLength= [0.2,inf]; % for final going to filter by the length of the bundle... 
        % this will automatically save short filo with embedded. 
        filterParams.saveFiloByLengthAndSig = []; 
        
        
    case 'ConnectToVeil_DensityOrient'; 
         filterParams.filoTypes = [0,1]; % 0 order attached to veil (no Branch), 1st order attached to a veil with a branch
         filterParams.filterByFit = 1; 
         filterParams.filterByBundleLength = [0.2,inf]; 
         
         filterParams.saveFiloByLengthAndSig = [5 inf; 50 100]; 
         % currently the first row is the min/max cut offs for the Ext Filo
         % the 2nd row is the upper and lower signal percentile 
         % Filopodia that fall into both these categories will be saved. 
         
         % for density you want to be slightly less rigid when filering out bad fits
         % bad fits can happen because of overlapping filo- want to include
         % these in the final calculation- do the same for orientation. 
   
    case 'Branch';
        % find maximum 
        filoTypes = [1 inf];  
        filterByFit = 1; 
end

%% START FILTERING
for iFrame = 1:length(analInfo)-1
    filoInfo = analInfo(iFrame).filoInfo;
    
    %% FILTER BY SIGMOIDAL FITTING
    if filterParams.filterByFit == 1;
        
        % find those data where the fit produces a NaN
        endpointCoordTest = vertcat(filoInfo(:).([ 'Ext_endpointCoordFitPix'])); % want this to filter out internal
        % put NaN for each filoInfo if the length metric is NaN.
        
        nanToRemove = isnan(endpointCoordTest);
        
        % find those data which did not pass the exit flag
        test = vertcat(filoInfo(:).Ext_exitFlag); % check exit flag - less than 1 means a majure failure
        % very infrequent though - can happen for a crossing check to see
        % what type of filopodia typical look like this.
        
        %logicals for filtering of filopodia
        toKeepBasedOnExit = (test>0 ); % good fits
        
        % eventually here want to create based on
    end % if filterByFit
    
    %% FILTER BY FILO TYPE
    
    type = vertcat(filoInfo(:).('type'));
    if filterParams.filoTypes(2) == inf
        % find the max branch type for that frame 
        filoTypesC = filterParams.filoTypes(1):max(type); 
    else 
        filoTypesC = filterParams.filoTypes; 
    end 
    % make cell array of logicals
    toKeepBasedOnTypeCell = arrayfun(@(x) type == filoTypesC(x) ,1:length(filoTypesC),'uniformoutput',0);
    
    % combine them
    toKeepBasedOnType = horzcat(toKeepBasedOnTypeCell{:});
    toKeepBasedOnType = sum(toKeepBasedOnType,2)~=0;
    %% FILTER BY TOTAL LENGTH 
    if ~isempty(filterParams.filterByBundleLength)
        minLength = filterParams.filterByBundleLength(1); 
        maxLength = filterParams.filterByBundleLength(2); 
        % get all lengths 
        lengthExt = vertcat(filoInfo(:).Ext_length); 
        lengthInt = vertcat(filoInfo(:).Int_length); 
        % change NaNs into zero for addition 
        lengthExt(isnan(lengthExt))= 0 ; 
        lengthInt(isnan(lengthInt)) = 0 ;
        lengthBundle = lengthExt + lengthInt; 
        % convert to um 
        lengthBundle = lengthBundle.* 0.216; % make input! 
        
        toKeepBasedOnLength = lengthBundle>minLength & lengthBundle<maxLength; % logical marking all the filo that 
        % make the bundle length criteria 
    end  % isempty 
%% Save Long Filopodia with Strong Signal If Desired
% note sometimes filopodia that cross will have a poor fit, these need to be 
% excluded from certain measurements (for instance length) and need
% to be included in others (for instance in a density metric). 

   if ~ isempty(filterParams.saveFiloByLengthAndSig); 
       s = filterParams.saveFiloByLengthAndSig;
       % get filopodia that meet length cut-off.. 
       lengthExt = vertcat(filoInfo(:).Ext_length)*.216; % convert 
       
       savePop1 = lengthExt>s(1,1) & lengthExt < s(2,1) ; 
       
       % get the full population 
       intensities = vertcat(filoInfo(:).Ext_IntensityNormToVeil); 
       intensitiesForPer = intensities(~isnan(intensities)); 
       cutoffMin = prctile(intensitiesForPer,s(2,1)); 
       cutoffMax = prctile(intensitiesForPer,s(2,2)); 
       
       
       savePop2 = intensities>cutoffMin & intensities<cutoffMax & ~isnan(intensities); 
       
       savePop = savePop1 & savePop2;
       % get filopoida that meet intensity cut-off (defined by percentile) 
     
   end 
    %% Make Final Filo Filter Set Based on All the Above Criteria 
    filoFilter = (toKeepBasedOnExit & toKeepBasedOnType & ~nanToRemove & toKeepBasedOnLength);
    if ~isempty(filterParams.saveFiloByLengthAndSig);  
        filoFilter = (filoFilter | savePop);  
    end
    
    % 
    filoFilterSet{iFrame} = filoFilter; 
     
    
end % for iFrame


end

