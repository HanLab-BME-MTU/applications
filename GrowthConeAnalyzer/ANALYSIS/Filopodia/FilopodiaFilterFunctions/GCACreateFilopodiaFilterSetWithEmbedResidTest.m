function [ filoFilterSet,filterParams] = GCACreateFilopodiaFilterSetWithEmbed(analInfo,filterType,varargin)
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
%       filoFilterSet:   rx2 cell of logical filters for the filopodia for each frame
%                        where r1 is the number of frames and column 1 is
%                        the traditional non-embedded filter and column2 is
%                        the embedded actin bundle filter. 
%
%
%
%% Check Input
ip = inputParser;
ip.CaseSensitive = false;
ip.addParameter('pixelSizeNm',216,@(x) isscalar(x));
ip.parse(varargin{:});
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
        filterParams.filterBasedOnGroupUpstream = 0; 
        filterParams.filterBasedOnBranchConnection =0 ; 
        filterParams.filterIntNoConnect=0;   
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
        filterParams.filterBasedOnGroupUpstream = 0; 
        filterParams.filterBasedOnBranchConnection =0 ; 
        filterParams.filterIntNoConnect=0;   
        
    case 'Branch2ndOrder_LengthInt';
        %
        filterParams.filoTypes = 2;
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.3, inf];
        filterParams.saveFiloByLengthAndSig = [];
        filterParams.embedFitCriteria = 95 ; % Not relavent
        filterParams.filoFitCriteria = 95;
        filterParams.filterBasedOnGroupUpstream = 0;  
        filterParams.filterBasedOnBranchConnection =0 ; 
        filterParams.filterIntNoConnect=0;   
        
    case 'Branch2ndOrder_Density_WithZeroOrder'; % WORKING
        filterParams.filoTypes = [0,1,2]; % 1st order attached to a veil with a branch, 2 branch
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.3,inf];
        %filterParams.saveFiloByLengthAndSig = [5 inf; 50 100];
        filterParams.saveFiloByLengthAndSig = [];
        filterParams.embedFitCriteria = 0 ; %
        filterParams.filoFitCriteria = 95;
        filterParams.filterBasedOnGroupUpstream = 0; 
        filterParams.filterBasedOnBranchConnection =0 ; 
        filterParams.filterIntNoConnect=0;   
        
    case 'Branch2ndOrder_Density_NoZero'; % WORKING
        filterParams.filoTypes = [1,2]; % 1st order attached to a veil with a branch, 2 branch
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.3,inf];
        filterParams.saveFiloByLengthAndSig = [5 inf; 50 100];
        filterParams.embedFitCriteria = 0 ; %
        filterParams.filoFitCriteria = 95;
        filterParams.filterBasedOnGroupUpstream = 0; 
        filterParams.filterBasedOnBranchConnection =0 ; 
        filterParams.filterIntNoConnect=0;   
        
    case 'Branch3rdOrder_LengthInt';
        filterParams.filoTypes = 3;
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.2, inf];
        filterParams.saveFiloByLengthAndSig = [];
        filterParams.embedFitCriteria = 0 ; %
        filterParams.filterBasedOnGroupUpstream = 0; 
        filterParams.filterBasedOnBranchConnection =0 ; 
        filterParams.filterIntNoConnect=0;   
        
        
    case 'Branch3rdOrder_DensityOrient';
        filterParams.filoTypes = 3; % 0 order attached to veil (no Branch), 1st order attached to a veil with a branch
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.2,inf];
        
        filterParams.saveFiloByLengthAndSig = [5 inf; 50 100];
        filterParams.embedFitCriteria = 0 ; %
        filterParams.filterBasedOnGroupUpstream = 0;
        filterParams.filterBasedOnBranchConnection =0 ; 
        
        filterParams.filterIntNoConnect=0;   
        
        
    case 'Branch4thOrder_LengthInt';
        filterParams.filoTypes = 4;
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.2, inf];
        filterParams.saveFiloByLengthAndSig = [];
        filterParams.filterBasedOnGroupUpstream = 0; 
        filterParams.filterBasedOnBranchConnection =0 ; 
        filterParams.filterIntNoConnect=0;   
        
        
    case 'Branch4thOrder_DensityOrient';
        filterParams.filoTypes = 3; % 0 order attached to veil (no Branch), 1st order attached to a veil with a branch
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.2,inf];
        
        filterParams.saveFiloByLengthAndSig = [5 inf; 50 100];
        filterParams.filterBasedOnGroupUpstream = 0; 
        filterParams.filterBasedOnBranchConnection =0 ; 
        filterParams.filterIntNoConnect=0;   
        
        
    case 'Validation';
        filterParams.filoTypes = [0 Inf]; % 0 order attached to veil (no Branch), 1st order attached to a veil with a branch
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.3,inf];
        
       % filterParams.saveFiloByLengthAndSig = [];
          filterParams.saveFiloByLengthAndSig = [1 inf; 5 100];
      
        filterParams.embedFitCriteria =95 ; % change threshold to be a percentile of the residuals 
        filterParams.filoFitCriteria = 95; 
        
        
        filterParams.filterBasedOnGroupUpstream = 0; 
        filterParams.filterBasedOnBranchConnection = 0 ; 
        
        filterParams.filterIntNoConnect=0; 
        
        % current problem is that need an option to maintain the flag on
        % the external but kick out the internal if the fit doesn't work
        % one option would be to kick out that info entirely (ie filter at
        % the level of the filoInfo fitting function (ie set these to NaNs)
        % maybe that is the better option just because this will save info
        % and if you do this here it will allow the user to make this more 
        % permissive or less.      
    case 'Validation2'  
        
        
              
   
        filterParams.filoTypes = [0 Inf]; % 0 order attached to veil (no Branch), 1st order attached to a veil with a branch
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.3,inf];
        
       filterParams.saveFiloByLengthAndSig = [];
          %filterParams.saveFiloByLengthAndSig = [5 inf; 50 100];
      
        filterParams.embedFitCriteria =95 ; % change threshold to be a percentile of the residuals 
        filterParams.filoFitCriteria = 95; 
        
        
        filterParams.filterBasedOnGroupUpstream = 1; 
        filterParams.filterBasedOnBranchConnection =1 ; 
        
        filterParams.filterIntNoConnect=0; 
        
        % current problem is that need an option to maintain the flag on
        % the external but kick out the internal if the fit doesn't work
        % one option would be to kick out that info entirely (ie filter at
        % the level of the filoInfo fitting function (ie set these to NaNs)
        % maybe that is the better option just because this will save info
        % and if you do this here it will allow the user to make this more 
        % permissive or less.       
    case 'Validation_NoEmbed'
        filterParams.filoTypes = [0 Inf]; % 0 order attached to veil (no Branch), 1st order attached to a veil with a branch
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.3,inf];
        
        filterParams.saveFiloByLengthAndSig = [];
        %filterParams.saveFiloByLengthAndSig = [5 inf; 50 100];
        filterParams.filoFitCriteria = 95;
        
        
        filterParams.filterBasedOnGroupUpstream = 0;
        filterParams.filterBasedOnBranchConnection = 0 ;
        
        
        
        
    case 'curvVsLength'
        filterParams.filoTypes = 0; % 0 order attached to veil (no Branch), 1st order attached to a veil with a branch
        filterParams.filterByFit = 1;
        filterParams.filterByBundleLength = [0.3,inf];
        
        % filterParams.saveFiloByLengthAndSig = [];
        filterParams.saveFiloByLengthAndSig = [5 inf; 50 100];
        
        filterParams.embedFitCriteria = 95 ; % change threshold to be a percentile of the residuals
        filterParams.filoFitCriteria = 95;
        filterParams.filterBasedOnGroupUpstream = 0; 
             filterParams.filterBasedOnBranchConnection =0 ; 
           filterParams.filterIntNoConnect=0;    
    case 'ConnectToVeil_LengthInt_Biosensors'
        filterParams.filoTypes = 0; % 0 order attached to veil (no Branch), 1st order attached to a veil with a branch
        filterParams.filterByFit = 1; % turn on filtering by filopodia fit.
        filterParams.filterByBundleLength= [0.3,inf]; % for final going to filter by the length of the bundle...
        % this will automatically save short filo with embedded.
        filterParams.saveFiloByLengthAndSig = [];
        filterParams.embedFitCriteria = 95 ; % change threshold to be a percentile of the residuals 
        filterParams.filoFitCriteria = 95; 
        filterParams.filterBasedOnGroupUpstream = 0; 
        filterParams.filterBasedOnBranchConnection =0 ; 
        filterParams.filterIntNoConnect=0;   
        
        % filter based on local threshold around filopodia in FRET/Donor channel 
        filterParams.sigmaNumerator = 3; % multiple times the standard deviation of the background estimation 
        % for the numerator of a FRET ratio - typically FRET channel 
        
        filterParams.sigmaDenominator = 1; % multiple time the standard deviation of the background estimation
        % for the denominator of the FRET ratio = typically the Donor
        % channel - as this value should decrease upon signal activity 
        % this should likely not be quite as stringent 
        
        filterParams.backSampleSize = 20; % the sample size for the background has to be at least this large  
        filterParams.excludeCrosses = false;  %need to implement this... 
        
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
%    % Test for empty filoInfo
%    
% 
%     end 
% end 
   
   
   
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
 

if isfield(filterParams,'embedFitCriteria') 
 
 
p1Int  = arrayfun(@(x) filoInfo(x).Int_params(1),1:length(filoInfo),'uniformoutput',0)'; % get the amplitude of the fit
resids = arrayfun(@(x) filoInfo(x).Int_resid,1:length(filoInfo),'uniformoutput',0)'; 

%filtInt = (p1Int<filterParams.embedFitCriteria(1) | isnan(p1Int)); 
% for a quick test filter out the amplitude of signal that is less than the
% 95th percentile of the residuals. 
filtInt  = cellfun(@(x,y) (y < prctile(x,filterParams.embedFitCriteria) | isnan(y)),resids,p1Int); 
else 
    filtInt = false(length(filoInfo),1); % no filtering based on embedded
end 
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
        if isfield(filterParams,'embedFitCriteria')
            lengthInt = vertcat(filoInfo(:).Int_length);
        else 
            lengthInt = zeros(length(filoInfo),1);
        end
        % change NaNs into zero for addition
        lengthExt(isnan(lengthExt))= 0 ;
        lengthInt(isnan(lengthInt)) = 0 ;
        lengthInt(filtInt) = 0; 
        lengthBundle = lengthExt + lengthInt;
        % convert to um
        lengthBundle = lengthBundle.* ip.Results.pixelSizeNm./1000; % most intuitive to set the length of the bundle in um
        
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
    
    %% Filter Based on Group Upstream
        toKeepBasedOnGroupUpstream = ones(length(filoInfo),1); 
    if filterParams.filterBasedOnGroupUpstream
        % get all the first order filopodia filtered (type 1)
        % if filtered propogate the filter to all downstream branches by getting
        % the idx of the group ID.
        types = vertcat(filoInfo(:).type);
        typeID = unique(types);
        typeID(typeID==0) = [];  % remove all nonbranch filo
        %         typeID(typeID~=0) =[];
        nTypes= length(typeID)-1; % don't need to look at the last type as there is
        % no branches attached to this.
  %      for i = 2:nTypes
   if ~isempty(typeID) 
       for i = 1:nTypes
           idxBranchStruct= vertcat(filoInfo(:).type) == typeID(i);
           
           idxBadFitWithAttach = ((~toKeepBasedOnExit | nanToRemove) & idxBranchStruct);
           
           % get the groupIDs of all the body attached filopodia
           groupIDsToFilter = vertcat(filoInfo(idxBadFitWithAttach).groupCount) ;
           groupIDsToFilter = unique(groupIDsToFilter);
           
           if ~isempty(groupIDsToFilter)
               % want to filter out any filo branches with a nType greater
               % than the low confidence branch previously filtered
               % so for each filoBranch detection of the current frame 
               % ask if it is member of the groupIDstoFilter and if 
               % it has a ID greater than the branch stem with attachment 
               % under question (don't want to filter high confidence
               % members of the group that are closer to the cell
               % veil/stem)
               toKeepBasedOnUpstreamNType(:,i) = arrayfun(@(x) ~(ismember(filoInfo(x).groupCount, groupIDsToFilter)...
                   & filoInfo(x).type > typeID(i)),1:length(filoInfo));
           else
               toKeepBasedOnUpstreamNType(:,i) = ones(length(filoInfo),1); % keep them all
           end
       end % i = 1:nTypes
        
        % if there is no flag to remove keep 
        toKeepBasedOnGroupUpstream = (sum(toKeepBasedOnUpstreamNType,2) == nTypes); 
   else % no branches to check 
          toKeepBasedOnGroupUpstream = ones(length(filoInfo),1); 
   end 
        clear toKeepBasedOnUpstreamNType
    end % filterBasedOnGroupUpstream
    %% Test for Connection after filter 
       idxKeepBasedOnBranchConnection = ones(length(filoInfo),1); 
    if filterParams.filterBasedOnBranchConnection
        % get the branch IDs 
             
         % for each connectivity index test if conXYCoords (3) is significantly > than 
         % length if it is get the corresponding conIdx
         
         % initiate logicalfilter vector: assume initially keep all filo
         idxKeepBasedOnBranchConnection = ones(length(filoInfo),1); 
         
         if isfield(filoInfo,'conXYCoords') % if no branches in the frame will not have this field. 
         for i = 1:length(filoInfo)
             if ~isempty(filoInfo(i).conXYCoords) 
                 % ask if the attachment distance is < the length metric 
                 % if it is keep if not it will be set to zero 
                 d = filoInfo(i).conXYCoords(:,3); 
                 all = filoInfo(i).Ext_params; 
                 conID = filoInfo(i).conIdx; 
                 typeC = filoInfo(i).type; 
                 if isnan(all)
                     l = 0; 
                 else 
                     l = filoInfo(i).Ext_params(2); 
                 end 
                 
                 idxRemovePerFilo = d>(l+1);
                 if sum(idxRemovePerFilo)~=0
                     temp = conID(idxRemovePerFilo); 
                     % get decendents 
                     groupC = filoInfo(conID(idxRemovePerFilo)).groupCount; 
                     groupIDs = vertcat(filoInfo(:).groupCount); 
                     typeAll = vertcat(filoInfo(:).type); 
                     
                     idxRemoveDownstream = find(arrayfun(@(x) ismember(filoInfo(x).groupCount,groupC)...
                         & filoInfo(x).type>typeC,1:length(filoInfo))); 
                     
                     toRemove{i} = idxRemoveDownstream'; 
                
                 end
                   else toRemove{i} =[]; 
             end 
           
             
         end
         idxNumToRemoveAll = vertcat(toRemove{:});
         if ~isempty(idxNumToRemoveAll)
             idxKeepBasedOnBranchConnection(idxNumToRemoveAll) = 0;
         end
         
         clear toRemove
         end % isfield(filoInfo,'conXYCoords')
    end % filterBasedOnBranchConnection
    %% Quick fix for embedded filo that hit a cell edge before being reconnected with partner. 
    idxKeepBasedOnIntNoConnect =  ones(length(filoInfo),1);
    if isfield(filterParams,'filterIntNoConnect')
        if filterParams.filterIntNoConnect
            % get the end points of the two coords
            % if not within 3 pixels remove
            %         idxKeepBasedOnIntNoConnect = ones(length(filoInfo),1);
            y = zeros(length(filoInfo(:)),1);
            for i = 1:length(filoInfo);
                if ~isnan(filoInfo(i).Int_coordsXY);
                    test1 = filoInfo(i).Int_coordsXY(1,:);
                    test2 = filoInfo(i).Ext_coordsXY(1,:);
                    if ~isnan(test1)
                        y(i) = sqrt((test1(1)-test2(1))^2 +(test1(2)-test2(2))^2);
                        
                    end
                end
            end
            
            idxKeepBasedOnIntNoConnect = y<6;
            %          idxKeepBasedOnIntNoConnect  = idxKeepBasedOnIntNoConnect';
        end % if filterParams.filterIntNoConnect
    end % if isfield
    %% Add Biosensor Filter 
    % Initiate
    toKeepBasedOnNumSigAtActinBundEnd = ones(length(filoInfo),1);
    toKeepBasedOnNumBackSampleSize = ones(length(filoInfo),1);
     if isfield(filterParams,'sigmaNumerator'); 
          % ok note there is a bit of a design problem as I had these as variable names 
         % when collecting - for now they should be the same though for all
         % so keep here inflexible. 
         
         % collect the background numerator values and estimate a local
         % threshold based on a simple 3*std of the background gaussian. 
         backValuesNum = arrayfun(@(x) filoInfo(x).FRET_Detect.BackgroundValues,1:length(filoInfo),'uniformoutput',0) ; 
         
         thresholds = cellfun(@(x) nanmean(x) + filterParams.sigmaNumerator*nanstd(x),backValuesNum); 
         
         % for now just filter by the estimation at the end of the actin
         % bundle - make sure to flip the dimensions as saved the indices
         % of the full actin bundle fit from outside veil (traditional filo) to inside
         % veil. ie The first point should be the lowest. 
         filoValuesNum =  arrayfun(@(x) flip(filoInfo(x).FRET_Detect.Ext_values),1:length(filoInfo),'uniformoutput',0); 
        
         % unfortunately have to do this via a formal loop as NaN as a
         % index will error
         for i = 1:length(filoInfo)
             idxC = filoInfo(i).indicesFitFullBundle(1); % get the discrete index for the 
             % end of the full actin bundle (ordered from external to
             % internal) 
             
             if ~isnan(idxC)
                 % get the intensity at the end of the actin bundle and
                 % test if it is higher than the threshold 
                 filoValuesNumAtPtC = filoValuesNum{i}(idxC);
                 toKeepBasedOnNumSigAtActinBundEnd(i) = filoValuesNumAtPtC > thresholds(i);
             else
                 toKeepBasedOnNumSigAtActinBundEnd(i) = 0;
             end
             
         end % for i = 1:length(filoInfo)
      
         toKeepBasedOnNumBackSampleSize = cellfun(@(x) length(x) > filterParams.backSampleSize,backValuesNum); 
          toKeepBasedOnNumBackSampleSize =   toKeepBasedOnNumBackSampleSize'; 
     end 
    
    
        %% Make Final Filo Filter Set Based on All the Above Criteria
        
        
        filoFilter = (toKeepBasedOnExit & toKeepBasedOnType & ~nanToRemove & toKeepBasedOnLength & toKeepBasedOnGroupUpstream & idxKeepBasedOnBranchConnection ... 
             &  toKeepBasedOnNumSigAtActinBundEnd  & toKeepBasedOnNumBackSampleSize);
        if ~isempty(filterParams.saveFiloByLengthAndSig);
            filoFilter = (filoFilter | savePop);
        end
        filtInt = ~filtInt & idxKeepBasedOnIntNoConnect; 
        %
        filoFilterSet{iFrame} = [filoFilter filtInt]; % originally 1 marked internal filter now 1 marks internal to consider.
        % as is and 0 marks filter.
    else
        filoFilterSet{iFrame} = [];
        
        display(['No FiloInfo for frame ' num2str(iFrame)]);
    end % for iFrame


end

