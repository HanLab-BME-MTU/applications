function [ groupData] = mitoticMeasureAndCollectMTLateralClassificationInfo(groupList,varargin)
% mitoticMeasureAndCollectMTLateralClassificationInfo: Measures the information
%   necessary to make the lateral classification for each MT track in a 
%   given subRegion as specified by groupList. Set up as a wrapper. If user
%   specifies it will collect this information in a groupData structure. 
%                                            
%                                            
% From plusTipTracker output : loads :
%                              /sub_X/meta/projData.mat
%
% From mitoticGetCellEdgeViaThresholding output: loads :
%                                                .tif files from 'masks':
%                                                (These are perFrame masks of the cell edge)
%                              
% INPUT:
% groupList: (REQUIRED) : Rx2 cell array
%                        holding directories of the subRegional MT tracks
%                        you would like to include is stored.
%
% minDispVect: (PARAM) : vector
%                        Default : 3 pixels
%                        A vector of min MT Track Lengths through which
%                        you would like to loop. 
%                        (currently in pixels)
%
% collectGroupData (PARAM) : logical
%                        Default : false 
%                        Collects the groupData file for further pooling
%                        and plotting of the 
%                        OrientVsDispVsDispMTVect_(ip.Results.minDispVect) 
%                        matrix
%
% groupDataFilename (PARM)  : character array or []
%                        Default : 'groupData'
%                        Provides a name for the groupData file. 
%                        Only applicable if collectGroupData is set to true
%                        
% outputDirectory: (PARAM) : character array
%                        Default: pwd (ie. current directory)
%                        The directory where the groupData.mat file will be
%                        saved if the collectGroupData option is turned on. 
%                        Only applicable if collectGroupData is set to true
%
% TSOverlays (PARAM)    : logical
%                         Default : true 
%                         Flag to make the troubleshooting (TS) overlays 
%                         that plot the classification measurements per track 
%                         extracted: note this will
%                         incease the time it take to run the results but 
%                         can be very helpful if you want to reset the
%                         classification cut-offs to target a specific biological 
%                         behavior.  
%% OUTPUT:
%  Primary Output:
%  For each project in the groupList creates and or adds a field to
%        projList{iProj,2} /meta/CorticalInfo/corticalData
%
%        added fields include
%              .OrientVsDispVsDispMTVect_(ip.Results.minDispVect)
%                       a rx3 double where r is the number of microtubules associated
%                       with a given subRegion and columns correspond to
%                       1) MT Orientation relative to cell edge (Orient)
%                       2) MT Cortical Distance (Disp)
%                       3) MT Track Length (DispMTVect)
%
%
%              .OrientCalcCompleteFlag_(minDispVect) 
%                       logical, set to 1 complete the measurements
%                       successfully 
%
%        Note: If the field already exists it will load these values but not
%        overwrite the measurements.
%
%  Optional Output:
%   If collectGroupData
%   groupData.mat file will be written to ip.Results.outputDirectory
%
%
%   If TSOverlays : 
%
%
%
% 
%% Check INPUT
ip = inputParser;

ip.CaseSensitive = false;

ip.addParameter('outputDirectory',pwd,@(x) ischar(x)); 

ip.addParameter('minDispVect',3,@(x) isnumeric(x));

ip.addParameter('collectGroupData',false); 

ip.addParameter('groupDataFilename','groupDataClassInfo');

ip.addParameter('TSOverlays',true); 


ip.parse(varargin{:});

%% Set up

%
if isempty(ip.Results.outputDirectory)
    saveDir = pwd;
else
    saveDir = ip.Results.outputDirectory;
end

if ~isdir(saveDir)
    mkdir(saveDir)
end
minDispVect = ip.Results.minDispVect;
btwGrpNames =  unique(groupList(:,1),'stable');
%% Start Loop
for i_minDisp = 1:length(minDispVect)
    minDispC = minDispVect(i_minDisp);
    count = 0;
    for iGroup = 1:length(btwGrpNames)
        projListIdxC = strcmpi(btwGrpNames{iGroup},groupList(:,1));
        projListC = groupList(projListIdxC,2);
        
        nProj = length(projListC(:,1));
        if ip.Results.collectGroupData
            groupData.names{iGroup} = btwGrpNames{iGroup};
            groupData.projList{iGroup} = projListC;
        end
        
        countBad = 1;
        for iProj = 1:nProj
            
            performOrientCalcs = 1; % assume in the beginning that you will perform the cals
            % get identifiers from filenames for .csv file
            it = 1;
            x = char(projListC(iProj));
            while it <5
                [x,id,num] = getFilenameBody(x);
                it = it +1;
            end
            y = char(projListC(iProj));
%             [~,sub,num2] = getFilenameBody(y);
%             
%             dataLabel{iProj+count} = [id num ' ' sub num2];
            s = load([projListC{iProj} filesep 'meta' filesep 'projData.mat']);
            projData =s.projData;
            
            % test quickly for number of growths
            %  if projData.stats.nGrowths == 0
            if isnan(projData.dwellAllTracks)
                performOrientCalcs = 0 ;
                display(['No Growth For Pole Cannot Perform Calcs for ' projListC{iProj} ]);
                if ip.Results.collectGroupData
                    groupData.noGrowth{iGroup}(countBad,1) = iProj;
                end
                countBad = countBad+1;
            end
            
            % if incomplete perform the orientation calculations and get a matrix of
            % the cortical classifiers to be used for that minDispVect (ie
            % orientation dispInRegin(same for each minDispVect) and the
            % actual MTDispVect used (as discrete)
            
            % load corticalDataMatifexists
            corticalDir = [projData.anDir filesep 'meta' filesep 'CorticalInfo'];
            if ~isdir(corticalDir)
                mkdir(corticalDir)
            end
            if exist([corticalDir filesep 'corticalData.mat'],'file')==2
                s2 = load([corticalDir filesep 'corticalData.mat']);
                corticalData = s2.corticalData;
                
                saveField = 0; % initiate a save field flag: default is to overwrite
                
                % check if there is a field that already exists
                if isfield(corticalData,['OrientVsDispVsDispMTVect_' num2str(minDispC)]);
                    % see if it was marked complete
                    if isfield(corticalData,['OrientCalcCompleteFlag_' num2str(minDispC)]);
                        saveField = corticalData.([ 'OrientCalcCompleteFlag_' num2str(minDispC)]); % if 1 it was complete flag to save (no overwrite).
                        
                    end
                    if saveField == 0 % if no flag to save field
                        warning(['Field for minDispVect' num2str(minDispC) 'Already Exists: Overwriting Field']); % give the user a warning
                        corticalData.(['OrientCalcCompleteFlag_' num2str(minDispC)]) = 0; % overwrite the fields
                        % re-set the current fields so can overwrite cleanly if necessary
                        corticalData.(['OrientVsDispVsDispMTVect_' num2str(minDispC)]) = [];
                    else
                        performOrientCalcs = 0; % field already calculated flag to skip calcs
                    end
                end % isfield(corticalData...)
            end % if exist(corticalData..)
            
            if performOrientCalcs == 1
                % add to corticalDataStructure.
                display(['Start Orient Calcs for' projListC{iProj}]);
                % note the wrapper set up is a bit redundant as it does
                % again load projData in mitoticLateralInfoClassMat...
                [OrientVsDispVsDispMTVect] = mitoticGetLateralInfoClassMat(projListC{iProj},... 
                    'minDispVect', minDispC,'TSOverlays',ip.Results.TSOverlays);
                % save and mark complete
                corticalData.(['OrientVsDispVsDispMTVect_' num2str(minDispC)]) = OrientVsDispVsDispMTVect; % put into a field might be a little better
                
                corticalData.(['OrientCalcCompleteFlag_' num2str(minDispC)]) = 1;
                
                save([corticalDir filesep 'corticalData.mat'],'corticalData');
            elseif countBad >1
                OrientVsDispVsDispMTVect = nan(1,3);
            else
                % load OrientVsDispVsDispMTVect
                OrientVsDispVsDispMTVect =   corticalData.(['OrientVsDispVsDispMTVect_' num2str(minDispC)]);
                
            end % if performOrientCalcs
            
            %% add to groupData set
            if ip.Results.collectGroupData
                
                groupData.(['OrientVsDispVsDispMTVect_' num2str(minDispC)]){iGroup}{iProj} = OrientVsDispVsDispMTVect;% name the field after the MTDispVect
                save([saveDir filesep ip.Results.groupDataFilename '.mat'],'groupData');
            end

            display(['Completed Running mDV Orientation Calcs for' projListC{iProj} ': mDV=' num2str(minDispC)  'Pixels' ]);
            clear corticalData
        end % iProj
        
        count = count + nProj;
        % clear the mats for next group
        clear orientAllProj dispAllProj dispVctMTAllProj
        display(['Completed Running mDV Orientation Calcs for' btwGrpNames{iGroup} ': mDV=' num2str(minDispC)  'Pixels' ]);
        
        % update GroupData
        if ip.Results.collectGroupData
            save([ip.Results.outputDirectory filesep ip.Results.groupDataFilename '.mat'],'groupData');
        end
    end % iGroup 

    display(['Completed Running all mDV Orientation Calcs for mDV=' num2str(minDispC)  'Pixels' ]);
end % i_minDisp

end
