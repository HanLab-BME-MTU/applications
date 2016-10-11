function [ filoFilterSet,filterParams] = GCACreateFilopodiaFilterSetWithEmbed(analInfo,filterType)
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
%      filterType:   PARAM:    char   Currently Either
%      'ConnectVeil_LengthInt',
%      'ConnectVeil_DensityOrient',
%      'Branch2ndOrder_LengthInt'
%      'Branch3rdOrder_LengthInt'
%      'Branch4thOrder_LengthInt'
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
        filterParams.filterByBundleLength= [0.3,inf]; % for final going to filter by the length of the bundle...
        % this will automatically save short filo with embedded.
        filterParams.saveFiloByLengthAndSig = [];
        filterParams.embedFitCriteria = 95 ; % change threshold to be a percentile of the residuals 
        filterParams.filoFitCriteria = 95; 
        
    case 'ConnectToVeil_DensityOrient';
        filterParams.filoTypes = [0,1]; % 0 order attached to veil (no Branch), 1st order attached to a veil with a branch
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.3,inf];
        
        filterParams.saveFiloByLengthAndSig = [5 inf; 50 100];
        
        % currently the first row is the min/max cut offs for the Ext Filo
        % the 2nd row is the upper and lower signal percentile
        % Filopodia that fall into both these categories will be saved.
        
        % for density you want to be slightly less rigid when filering out bad fits
        % bad fits can happen because of overlapping filo- want to include
        % these in the final calculation- do the same for orientation.
        
        filterParams.embedFitCriteria = 95 ; % change threshold to be a percentile of the residuals 
        filterParams.filoFitCriteria = 95; 
        
        
    case 'Branch2ndOrder_LengthInt';
        %
        filterParams.filoTypes = 2;
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.3, inf];
        filterParams.saveFiloByLengthAndSig = [];
        filterParams.embedFitCriteria = 95 ; % Not relavent
        filterParams.filoFitCriteria = 95;
         
        
    case 'Branch2ndOrder_Density_WithZeroOrder'; % WORKING
        filterParams.filoTypes = [0,1,2]; % 1st order attached to a veil with a branch, 2 branch
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.3,inf];
        %filterParams.saveFiloByLengthAndSig = [5 inf; 50 100];
        filterParams.saveFiloByLengthAndSig = [];
        filterParams.embedFitCriteria = 0 ; %
        filterParams.filoFitCriteria = 95;
        
        
    case 'Branch2ndOrder_Density_NoZero'; % WORKING
        filterParams.filoTypes = [1,2]; % 1st order attached to a veil with a branch, 2 branch
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.3,inf];
        filterParams.saveFiloByLengthAndSig = [5 inf; 50 100];
        filterParams.embedFitCriteria = 0 ; %
        filterParams.filoFitCriteria = 95;
        
    case 'Branch3rdOrder_LengthInt';
        filterParams.filoTypes = 3;
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.2, inf];
        filterParams.saveFiloByLengthAndSig = [];
        filterParams.embedFitCriteria = 0 ; %
        
    case 'Branch3rdOrder_DensityOrient';
        filterParams.filoTypes = 3; % 0 order attached to veil (no Branch), 1st order attached to a veil with a branch
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.2,inf];
        
        filterParams.saveFiloByLengthAndSig = [5 inf; 50 100];
        filterParams.embedFitCriteria = 0 ; %
        
        
        
    case 'Branch4thOrder_LengthInt';
        filterParams.filoTypes = 4;
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.2, inf];
        filterParams.saveFiloByLengthAndSig = [];
        
    case 'Branch4thOrder_DensityOrient';
        filterParams.filoTypes = 3; % 0 order attached to veil (no Branch), 1st order attached to a veil with a branch
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.2,inf];
        
        filterParams.saveFiloByLengthAndSig = [5 inf; 50 100];
        
    case 'Validation';
        filterParams.filoTypes = [0 Inf]; % 0 order attached to veil (no Branch), 1st order attached to a veil with a branch
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.3,inf];
        
       % filterParams.saveFiloByLengthAndSig = [];
          filterParams.saveFiloByLengthAndSig = [5 inf; 50 100];
      
        filterParams.embedFitCriteria = 95 ; % change threshold to be a percentile of the residuals 
        filterParams.filoFitCriteria = 95; 
        
        
        
        % current problem is that need an option to maintain the flag on
        % the external but kick out the internal if the fit doesn't work
        % one option would be to kick out that info entirely (ie filter at
        % the level of the filoInfo fitting function (ie set these to NaNs)
        % maybe that is the better option just because this will save info
        % and if you do this here it will allow the user to make this more 
        % permissive or less.       
    case 'curvVsLength'
        filterParams.filoTypes = 0; % 0 order attached to veil (no Branch), 1st order attached to a veil with a branch
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.3,inf];
        
        % filterParams.saveFiloByLengthAndSig = [];
        filterParams.saveFiloByLengthAndSig = [5 inf; 50 100];
        
        filterParams.embedFitCriteria = 95 ; % change threshold to be a percentile of the residuals
        filterParams.filoFitCriteria = 95;
        
        
        
end

%% START FILTERING
if length(analInfo) == 1 
    endFrame = 1; 
else 
    endFrame = length(analInfo)-1; 
end 

for iFrame = 1:endFrame
    filoInfo = analInfo(iFrame).filoInfo;
    
    if ~isempty(filoInfo); 
   % Test for empty filoInfo
   

   
   
   
   
%     noFiloInfoFlag =  arrayfun(@(x) isempty(analInfo(x).filoInfo),1:length(analInfo));
%     if sum(noFiloInfoFlag) ~=0
%         problems = find(noFiloInfoFlag);
%         arrayfun(@(x) display(['No filoInfo for Frame ' num2str(problems(x))]),1:length(problems));
%     end
    
%     for iFrame = 1:numel(analInfo) -1
%         filoInfo = analInfo(iFrame).filoInfo;
%         % arrayfun(@(x)
%         % test if the field associated with endpoint coordinate is empty
%         % 1 if missing fits 0 if not
%         missingFits(iFrame,1)= sum(arrayfun(@(x) isempty(x.Ext_endpointCoordFitPix),filoInfo));
%     end
 %% filter by embedded first 
p1Int  = arrayfun(@(x) filoInfo(x).Int_params(1),1:length(filoInfo),'uniformoutput',0)'; % get the amplitude of the fit
resids = arrayfun(@(x) filoInfo(x).Int_resid,1:length(filoInfo),'uniformoutput',0)'; 

%filtInt = (p1Int<filterParams.embedFitCriteria(1) | isnan(p1Int)); 
% for a quick test filter out the amplitude of signal that is less than the
% 95th percentile of the residuals. 
 filtInt  = cellfun(@(x,y) (y < prctile(x,filterParams.embedFitCriteria) | isnan(y)),resids,p1Int); 

    %% FILTER BY SIGMOIDAL FITTING
    if filterParams.filterByFit == 1;
        
        % find those data where the fit produces a NaN
        endpointCoordTest = vertcat(filoInfo(:).([ 'Ext_endpointCoordFitPix'])); % want this to filter out internal
        % put NaN for each filoInfo if the length metric is NaN.
        
        nanToRemove = isnan(endpointCoordTest);
        
        % find those data which did not pass the exit flag
        test = vertcat(filoInfo(:).Ext_exitFlag) ; % check exit flag - less than 1 means a majure failure
        % very infrequent though - can happen for a crossing check to see
        % what type of filopodia typical look like this.
        
        %logicals for filtering of filopodia
        toKeepBasedOnExit = (test>0 ); % good fits
        
        %filterParams.filoFitCriteria
        p1Ext  = arrayfun(@(x) filoInfo(x).Ext_params(1),1:length(filoInfo),'uniformoutput',0)'; % get the amplitude of the fit
        residsExt = arrayfun(@(x) filoInfo(x).Ext_resid,1:length(filoInfo),'uniformoutput',0)'; 

        %filtInt = (p1Int<filterParams.embedFitCriteria(1) | isnan(p1Int)); 
        % for a quick test filter out the amplitude of signal that is less than the
        % 95th percentile of the residuals. 
        filtExt  = cellfun(@(x,y) (y < prctile(x,filterParams.filoFitCriteria) | isnan(y)),residsExt,p1Ext); 
        toKeepBasedOnExit =  toKeepBasedOnExit & ~filtExt; 
        % eventually here want to create based on
    end % if filterByFit
    
    %% FILTER BY FILO TYPE
    
    type = vertcat(filoInfo(:).('type'));
    if length(filterParams.filoTypes)>1 && filterParams.filoTypes(2) == inf
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
        lengthInt(filtInt) = 0; 
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
        
        savePop1 = lengthExt>s(1,1) & lengthExt < s(2,1) & toKeepBasedOnType ; % fixed 201508
        
        % get the full population
        intensities = vertcat(filoInfo(:).Ext_IntensityNormToVeil);
        intensitiesForPer = intensities(~isnan(intensities));
        cutoffMin = prctile(intensitiesForPer,s(2,1));
        cutoffMax = prctile(intensitiesForPer,s(2,2));
        
        
        savePop2 = intensities>cutoffMin & intensities<cutoffMax & ~isnan(intensities) & toKeepBasedOnType; %fixed 201508
        
        savePop = savePop1 & savePop2;
        % get filopoida that meet intensity cut-off (defined by percentile)
        
    end
    %% Make Final Filo Filter Set Based on All the Above Criteria
    filoFilter = (toKeepBasedOnExit & toKeepBasedOnType & ~nanToRemove & toKeepBasedOnLength);
    if ~isempty(filterParams.saveFiloByLengthAndSig);
        filoFilter = (filoFilter | savePop);
    end
    
    %
    filoFilterSet{iFrame} = [filoFilter ~filtInt]; % originally 1 marked internal filter now 1 marks internal to consider. 
    % as is and 0 marks filter. 
    else 
        filoFilterSet{iFrame} = []; 
    
        display(['No FiloInfo for frame ' num2str(iFrame)]);
end % for iFrame


end

