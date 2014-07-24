function [ micropatternOutput ] = micropatternBtwGrpAnalysis(micropatternOutput,saveDir,makeHist)
%perfrom the between group comparisons for each mask type and extract type
%and save.

% Adds extra fields to micropatternOutput 
% micropatternOutput{iMask}.groupData{iExtract}{iRegion}{iWindow} if
% subregions 
% or 
% micropatternOutput{iMask}.groupData{iExtract}{iWindow} if no subregions; 

% order for subRegions is 1 = NonAd, 2 = Ad, 3 = Ad_Corn 
if nargin<1 || isempty(micropatternOutput)
    
    [file ,Path] =  uigetfile(pwd, 'Please Select MicropatternOutput File to Analysize')  ;
    s = load([Path filesep file]);
    micropatternOutput = s.micropatternOutput;
    
end

if nargin<2 || isempty(saveDir)
    saveDir =  uigetdir(pwd,'Please Select a Directory to Store The Analysis') ;
end

if nargin<3 || isempty(makeHist)
    makeHist = 0 ;
end 

%% Loop over Mask Type %%

for iMask = 1:numel(micropatternOutput)
    % load mask params
    
    maskParamsCurrent = micropatternOutput{iMask}.maskParams;
    groupLists = micropatternOutput{iMask}.groupList; % change in prepare
    winNum = maskParamsCurrent.numWindows;
    winSize = maskParamsCurrent.windowSize;
   
    % loop over the different extract types
    
    for iExtract = 1:numel(maskParamsCurrent.extractType)
        
        %% make general folders
        % get mask Identifies for folder names
        
        s1 = num2str(maskParamsCurrent.numWindows);
        s2 = num2str(maskParamsCurrent.windowSize);
        if maskParamsCurrent.subRegions == 1
            s3 = 'subRegions';
        else
            s3 = '';
        end ;
        
        maskDescrip = ['SUBROIS_' s1 '_' s2 '_umWindows_' char(s3) '_'];
        
        extractName = char(maskParamsCurrent.extractType{iExtract});
        
        
        folderMaskDescrip = [saveDir filesep maskDescrip] ;
        folderMaskDescripExtract = [folderMaskDescrip filesep extractName];
        
        if ~isdir([folderMaskDescripExtract filesep 'SubRegions'])
            mkdir([folderMaskDescripExtract filesep 'SubRegions'])
        end
        
        micropatternOutput{iMask}.primaryStatDir = folderMaskDescripExtract;
        
        %% make subregional folders and perform analysis
        if maskParamsCurrent.subRegions == 1
            
            % subregion folders
            %         if maskParamsCurrent.subRegions == 1
            %            fieldnames(groupData.pooledStats );
            %           for iField = 1:length(fieldnames)
            %
            
            
            region{1,1} = 'NonAd_';
            region{2,1} = 'Ad_';
            region{3,1} = 'AdCorn_';
            % make a list of all groupLists we want to analyze
            regLong{1,1} = 'NonAdhesion';
            regLong{2,1} = 'Adhesion';
            regLong{3,1} = 'AdhesionCorner';
            
            
            
            for iReg = 1:3 % for H-pattern will always have 3 regions
                for iWind = 1:winNum+1; % windows
                    if iWind == winNum+1
                        add = 'GreaterThan_';
                        dist  = (iWind-1)*winSize;
                    else
                        add = '';
                        dist = iWind*winSize;
                    end
                    
                    %     dir  = [folderMaskDescripExtract filesep 'SubRegions' filesep regLong{iReg,1}];
                    %
                    %
                    %     if ~isdir(dir)
                    %     mkdir(dir)
                    %     end
                    %
                    dir2 = [folderMaskDescripExtract filesep 'SubRegions' filesep regLong{iReg,1} filesep add num2str(dist) 'uM'] ;
                    if ~isdir(dir2)
                        mkdir(dir2)
                    end
                    
                    % load groupList
                    groupList = groupLists{iExtract}{iReg}{iWind};
                    save([dir2 filesep 'groupList'],'groupList'); 
                    
                  
                    
                    % make new directories
                    if ~isdir([dir2 filesep 'perCell'])
                        mkdir([dir2 filesep 'perCell'])
                    end
                    
                    if ~isdir([dir2 filesep 'pooled'])
                        mkdir([dir2 filesep 'pooled'])
                    end
                    
                    % perform the analysis
                    
                    % collect groupData
                    groupDataCurrent = plusTipExtractGroupData_workingMB(groupList);
                    save([dir2 filesep 'groupData'],'groupDataCurrent');
                    micropatternOutput{iMask}.groupData{iExtract}{iReg}{iWind} = groupDataCurrent;
                    micropatternOutput{iMask}.statDirInd{iExtract}{iReg}{iWind} = dir2; 
                    
                    % perform analysis 
                    plusTipGetHits(groupDataCurrent, [dir2 filesep 'perCell'], 0.05, 1,20);
                    plusTipPoolGroupData(groupDataCurrent,[dir2 filesep 'pooled'],1,makeHist);
                    plusTipTestDistrib(groupDataCurrent,[dir2 filesep 'pooled'],0.05,20,21); %
                    
                end % iWind
            end % iRegion
        else  % make folders for nonregional subtype (slightly redundant)
            
            numRois = maskParamsCurrent.numWindows +1;
            
            for iSubRoi = 1:numRois
                
                dir2 = [folderMaskDescripExtract filesep 'SubRegions' filesep 'sub_' num2str(iSubRoi)];
                
                if ~isdir(dir2)
                    mkdir(dir2)
                end
                
                groupList = groupLists{iExtract}{iSubRoi};
                
                save([dir2 filesep 'groupList'], 'groupList') % save a copy of the groupList in folder
               
                
                
                % make folders
                
                if ~isdir([dir2 filesep 'perCell'])
                    mkdir([dir2 filesep 'perCell'])
                end
                
                if ~isdir([dir2 filesep 'pooled'])
                    mkdir([dir2 filesep 'pooled'] )
                end
                
                 
                % collect and save groupData 
                groupDataCurrent = plusTipExtractGroupData_workingMB(groupList); % extact groupData
                save([dir2 filesep 'groupData'],'groupDataCurrent');
                micropatternOutput{iMask}.groupData{iExtract}{iSubRoi} = groupDataCurrent; % saveGroupData in Micropattern Output
                micropatternOutput{iMask}.statDir{iExtract}{iSubRoi} = dir2;
                
                % perform analysis 
                plusTipGetHits(groupDataCurrent, [dir2 filesep 'perCell'], 0.05, 1,20);
                plusTipPoolGroupData(groupDataCurrent,[dir2 filesep 'pooled'],1,makeHist);
                plusTipTestDistrib(groupDataCurrent,[dir2 filesep 'pooled'],0.05,20,21); %
                
                
                
                
            end   % iSubRois
            
            
            
        end % if maskParamsCurrent.subRegions == 1
        
        
        
        
        
        
        
        
        
        
        
    end % iExtract
    
end % iMask 

save([saveDir filesep 'micropatternOutput.mat'],'micropatternOutput'); % resave new file


end







%



