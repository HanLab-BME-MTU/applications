function [ groupData,noGrowthList] = newCorticalClassworking(groupList,saveDir,minDispVect,dayName)
%Function loops through the MT displacement vectors used for the
%orientaiton calcultion 
% INPUT: groupList of names and project paths 
%        saveDir: where the groupData will be saved
%        minDispVect: a vector of minDisp values one would like to loop
%        through : note for now this is designed to be a one time deal. 
%        dayName: ID for the groupData
% 
btwGrpNames =  unique(groupList(:,1),'stable');

for i_minDisp = 1:length(minDispVect)
    minDispC = minDispVect(i_minDisp);
    count = 0;
    for iGroup = 1:length(btwGrpNames)
        projListIdxC = strcmpi(btwGrpNames{iGroup},groupList(:,1));
        projListC = groupList(projListIdxC,2);
       
        
      
        
        nProj = length(projListC(:,1));
        groupData.names{iGroup} = btwGrpNames{iGroup};
        groupData.projList{iGroup} = projListC;
         
     
        
        
        % initialize a matrix where rows = MT track and columns = project number
        % (overinitialize rows with NaN  so can concatenate easily)
        orientAllProj = nan(300,length(projListC));
        dispAllProj = nan(300,length(projListC));
        dispVectMTAllProj = nan(300,length(projListC)); 
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
            [~,sub,num2] = getFilenameBody(y);
            
            dataLabel{iProj+count} = [id num ' ' sub num2];
            s = load([projListC{iProj} filesep 'meta' filesep 'projData.mat']);
            projData =s.projData;
            
            % test quickly for number of growths 
          %  if projData.stats.nGrowths == 0
          if isnan(projData.dwellAllTracks)
          performOrientCalcs = 0 ;
                display(['No Growth For Pole Cannot Perform Calcs for ' projListC{iProj} ]); 
                groupData.noGrowth{iGroup}(countBad,1) = iProj; 
                
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
               if isfield(corticalData,'OrientVsDispVsDispMTVect'); 
                   % remove old fields from previous iterations
                   corticalData = rmfield(corticalData,'OrientVsDispVsDispMTVect'); 
                   corticalData = rmfield(corticalData,'MinDisp'); 
               end 
              
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
              
                
            %   corticalData.MinDisp(i_minDisp) = 0; 
               
            end % if exist(corticalData..)
           
            if performOrientCalcs == 1
                % add to corticalDataStructure.
                display(['Start Orient Calcs for' projListC{iProj}]);
                [OrientVsDispVsDispMTVect] = plotIncidenceVsDispVsKappaNew(projListC{iProj},minDispC);
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
            groupData.(['OrientVsDispVsDispMTVect_' num2str(minDispC)]){iGroup}{iProj} = OrientVsDispVsDispMTVect;% name the field after the MTDispVect
            save([saveDir filesep 'groupData.mat'],'groupData');
            orientVectC = OrientVsDispVsDispMTVect(:,1);
            dispVectC = OrientVsDispVsDispMTVect(:,2);
            dispVectMTC = OrientVsDispVsDispMTVect(:,3); 
            %         % save individual project values in larger matrix 
            orientAllProj(1:length(orientVectC),iProj)= orientVectC;
            dispAllProj(1:length(dispVectC),iProj) = dispVectC;
            dispVectMTAllProj(1:length(dispVectMTC),iProj) = dispVectMTC; 
            
             
             display(['Completed Running mDV Orientation Calcs for' projListC{iProj} ': mDV=' num2str(minDispC)  'Pixels' ]);    
           clear corticalData
        end % iProj
        
        
        % save to groupData
        groupData.(['orientMatAllProj_' num2str(minDispC)]){iGroup} = orientAllProj;
        groupData.(['dispMatAllProj_' num2str(minDispC)]){iGroup} = dispAllProj;
        groupData.(['dispMTVectMatAllProj_' num2str(minDispC)]){iGroup} = dispVectMTAllProj; 
        count = count + nProj;
        % clear the mats for next group
        clear orientAllProj dispAllProj dispVctMTAllProj
    display(['Completed Running mDV Orientation Calcs for' btwGrpNames{iGroup} ': mDV=' num2str(minDispC)  'Pixels' ]); 
    
    % update GroupData
    save([saveDir filesep 'groupData.mat'],'groupData');
    end 
    % now collect all the columns of classifier data for each project over
    % all groups for export to a .csv file. 
    orientMatAll = horzcat(groupData.(['orientMatAllProj_' num2str(minDispC)]){:});
    dispMatAll = horzcat(groupData.(['dispMatAllProj_' num2str(minDispC)]){:});
    dispMTVectMatAll = horzcat(groupData.(['dispMTVectMatAllProj_' num2str(minDispC)]){:}); 
    
    toExport = mat2dataset(orientMatAll,'VarNames',dataLabel);
    %save([saveDir filesep [dayName 'velocityMat'],'velMat']);
    export(toExport,'file',[saveDir filesep dayName (num2str(minDispC)) '_cortClassOrient_degrees.csv']);
    toExport2 = mat2dataset(dispMatAll,'VarNames',dataLabel);
    export(toExport2, 'file',[saveDir filesep dayName  num2str(minDispC) '_cortClassDisp_um.csv']) ; % should be the same for all minDispC just a good sanity check. 
    toExport3 = mat2dataset(dispMTVectMatAll,'VarNames',dataLabel); 
    export(toExport3,'file',[saveDir filesep dayName num2str(minDispC) '_dispMTVect_pixels.csv']); 
    clear toExport toExport2 toExport3 orientMatAll dispMatAll dispMTMatAll
    
     display(['Completed Running all mDV Orientation Calcs for mDV=' num2str(minDispC)  'Pixels' ]); 
end % i_minDisp

end
