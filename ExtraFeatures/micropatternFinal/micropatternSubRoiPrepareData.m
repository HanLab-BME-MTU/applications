function [ micropatternOutput] = micropatternSubRoiPrepareData(groupList,maskParams,saveDir,varargin)
% micropatternSubRoiPrepareData : Function automatically partitions the
% HPattern Data Into Different SubRegions (Ad,NonAd,AdCorn)
%
%
% This function makes the respective subrois given the maskParams input for
% all projects in the groupList and does 1 or more subtrack extractions
% based for each type of subRoiMasking as specificid by
% maskParams.extractType
% It also automatically makes the subRoiGroupLists required for subsequent
% analysis
%
% INPUT:
% groupList: cell array of names and project directories
% NOTE please do not include numbers in the GroupName as we cannot yet
% accomodate this when making the groupLists
%
% maskParams: cell array of structures containing the mask parameters
%             for each mask with fields
%             maskParams{iMask}
%             .numWindows (number of regularly sized windows)
%             .windowSize (size of window)
%             .subRegions  1 make the extra regions, 0 do not
%             .extractType{j}  (can enter multiple extract types and will
%             loop through each for each iMask.
%             Choices for extraction: 'subTrackEnd','nonDirectional','nucleation'
%             If empty will load example input.
%
% saveDir : final directory to save the micropattern output file: if
%           saveDir = [] will ask to choose directory.
%
%  micropattern output: structure of cell arrays such that
%  micropatternOutput{iMask} with respective fields groupList and
%  corresponding maskParams
%  each groupList is a cell array of the form
%  groupList{iExtractType}{iSubRegion (if Applicable example non-ad, ad,
%  and corn)}{iWind (1 is closest to cell edge)}
%
% OUTPUT: micropatternOutput numelements = the cell number of masks in an
% analysis bundle
% structure containing two fields
% micropatternOutput.groupListsForAnalysis.btw = cell array containing all groupLists for comparing the different regions of the cell between different groups
% {number of
% micropatternOutput.maskParams = saves the mask params so as identifies
% for the group Lists
%%
% Input check
ip = inputParser;
ip.addRequired('groupList',@(x)iscell(x) || isempty(x));
ip.addRequired('maskParams',@(x)iscell(x) || isempty(x));%
ip.addRequired('saveDir', @(x)ischar(x) || isempty(x)); % eventually should make optional
ip.addOptional('append', 0 ,@isnumeric);

ip.parse(groupList,maskParams,saveDir,varargin{:});

append = ip.Results.append;

% nice to be able to append a micropattern array if want to add more
% micropatterns
if append ==1
    [filenameAppend, path] = uigetfile(pwd, 'Please Choose a Micropattern File to Append');
    s.micropatternOutput = load([path filesep filenameAppend]);
    micropatternOutput{:} = s.micropatternOutput{:};
else % do not append
end

%% % example maskParamInput (default)
%'
if isempty(groupList) 
 [file, path] = uigetfile(pwd,'Please Select the groupList.mat file you would like to use'); 
 s = load([path filesep file]); 
 groupList = s.groupList; 
end 

if isempty(maskParams)
  [file,path]  = uigetfile(pwd,'Please Select the maskParams.mat file you would like to use');
    s = load([path filesep file]); 
    maskParams = s.maskParams;
end 



%     maskParams{1}.numWindows = 1;
%     maskParams{1}.windowSize = 1;
%     maskParams{1}.subRegions = 1;
%     maskParams{1}.extractType{1} = {'subTrackEnd'} ;
%     %maskParams{1}.extractType{2} = {'subTrackStart'} ;
%     %maskParams{1}.extractType{3} = {'nucleation'};


if isempty(saveDir)
    saveDir =  uigetdir(pwd,'Please Select A Directory to Save the Output');
end
%% Collect Group Names

projGroupName = groupList(:,1);
%projGroupDir=cellfun(@(x) formatPath(x),groupList(:,2),'UniformOutput',0);
projGroupDir =groupList(:,2); % try not formating path for now.

% fix the names if there are spaces or hyphens and append prefix 'grp'
projGroupName=cellfun(@(x) ['grp_' regexprep(x,'[ -]','_')],...
    projGroupName,'uniformoutput',0);

% count unique groups and keep them in order of the original list
[btwGrpNames,m] = unique(projGroupName);
[~,idx]=sort(m);
btwGrpNames=btwGrpNames(idx);



%% Loop Through all Masks

for iMask =  1:numel(maskParams)
    
    clear groupLists % should reinitiate below but just in case.
    
    
    % get current mask parameters
    maskParamsCurrent = maskParams{iMask};
    
    
    
    
    numExtractTypes = numel(maskParamsCurrent.extractType);
    
    % reinitiate groupLists
    if maskParamsCurrent.subRegions == 1
        numRois = (maskParamsCurrent.numWindows+1)*3; % assume H pattern
        groupLists = cell(1,numExtractTypes);
        groupListsSPreCluster = cell(1,numExtractTypes);
    else
        numRois = maskParamsCurrent.numWindows+1;
        groupLists = cell(1,numExtractTypes);
        groupListsNSPreCluster = cell(1,numExtractTypes);
    end
    
    
    
    
    
    %     %upThree = getFilenameBody(upTwo);
    %     s1 = num2str(maskParamsCurrent.numWindows);
    %     s2 = num2str(maskParamsCurrent.windowSize);
    %     if maskParamsCurrent.subRegions == 1
    %         s3 = 'subRegions';
    %     else
    %         s3 = '';
    %     end ;
    %
    %     maskDescrip = ['SUBROIS_' s1 '_' s2 '_umWindows_' char(s3) '_'];
    %
    
    
    for iGroup = 1:length(btwGrpNames)
        projIndx=find(strcmp(btwGrpNames(iGroup),projGroupName));
        groupName = btwGrpNames(iGroup);
        % indices of projects in iGroup
        badMaskCount = 1; % reset badMask count for each group
        badMaskIdx = [];
        for iIndx = 1:length(projIndx)
            flag2stop = 0 ; %reset badMask Flag
            
            iProj = projIndx(iIndx);
            currentProjDir = char(projGroupDir(iProj));
            
            
            
            % Extract all masks specified by user for a given project
            [roiSet roiYX] = subRoiMakeMasksHPattern(currentProjDir,maskParamsCurrent );
            
            if iIndx ==1
                
                exampleSubRoiSet= roiSet;
                exampleRoiYXWholeCell = roiYX;
                
            end
            
            % loop over each type of subtrack extraction process
            % loop at this step as the same mask requires different
            % types of extractions and don't want to make the mask mult
            % times
            % also do not have to save and then reload the subRoi masks.
            % however need a pointer to the directory this is
            
            
            for iExtract = 1:numel(maskParamsCurrent.extractType)
                if flag2stop == 1
                    break
                end
                toExtract = maskParamsCurrent.extractType{iExtract};
                
                % perform each type of extraction
                
                % save the folder for each extract type
                
                
                % this will make the folders for the project
                [groupListIndProj subRoiDirs] = makeMicropatternFolders(currentProjDir,iIndx,roiSet,maskParamsCurrent,groupName,toExtract);  %%%
                
                if isempty(groupListIndProj) %  check for a bad mask will produce this as empty
                    
                    maskParamsCurrent.badMasks{badMaskCount} = currentProjDir;
                    badMaskIdx(badMaskCount) = iIndx; % get index of project %
                    % don't know how many there will be
                    badMaskCount = badMaskCount + 1; % waste of computational time but ok for a hack job
                    flag2stop= 1 ;
                    continue  % skip the rest
                    
                end
                
                
                if maskParamsCurrent.subRegions == 0
                    
                    groupListsNSPreCluster{iExtract}{iGroup}{iProj} = groupListIndProj;
                else % get lists for subregions below (for historical reasons easier) probably not the ideal way to have implemented
                    
                end
                
                % Extract data from each subregion
                for iSubDir = 1:numel(subRoiDirs)
                    subRoiDir = char(subRoiDirs(iSubDir)) ;
                    
                    plusTipSubRoiExtractTracksUpdate(subRoiDir,'fraction',0.10,char(toExtract)); % NEED TO fix not ideal input
                end
            end %
            
        end % for iProj
        
        
        
        
        % for each extract type
        for iExtract = 1:numel(maskParamsCurrent.extractType)
            if maskParamsCurrent.subRegions == 1
                % make withingroup proj lists
                projList = projGroupDir(projIndx);
                if ~isempty(badMaskIdx)
                    projList(badMaskIdx,:) = []; % don't include projects with badMasks
                end
                
                extractType = char(maskParamsCurrent.extractType{iExtract});
                if ~isempty(projList)
                    [listsForCurrentGroup] = makeGroupListsMicropatternWithinGroup(projList,maskParams{iMask},extractType) ;
                    groupListsSPreCluster{iExtract}{iGroup} = listsForCurrentGroup;
                else
                    groupListsSPreCluster{iExtract}{iGroup} = [];
                end
            else
            end
        end
    end % for each Group
    
    %% Make Between Group comparisons Lists
    
    %     folderMask = [btwGrpListSaveDir filesep maskDescrip];
    %
    %     if ~isdir(folderMask)
    %         mkdir(folderMask);
    %     end
    
    
    
    
    
    if maskParamsCurrent.subRegions ~=1
        
        for iExtract = 1:numel(maskParamsCurrent.extractType)
            %extractType = char(maskParamsCurrent.extractType{iExtract});
            
            %             folderMaskExtract = [folderMask filesep extractType];
            
            
            %             if ~isdir(folderMaskExtract)
            %                 mkdir(folderMaskExtract)
            %             end
            
            
            % cluster groupList
            btwGrpComp = cell(length(btwGrpNames),1);
            for iGroup = 1:length(btwGrpNames)
                btwGrpComp{iGroup,1} = vertcat(groupListsNSPreCluster{iExtract}{iGroup}{:});
            end
            
            % Get all group Information for that Mask/Extract Type
            combine = vertcat(btwGrpComp{:,1});
            
            % make individual subRoi groupLists % note if have more than 10
            % subrois this technique is not useful..however this will not be
            % the case for us
            for i = 1:numRois
                idxRegion = cellfun(@(x) ~isempty(strfind(x,num2str(i))),combine(:,1),'UniformOutput',0);
                groupLists{iExtract}{i} = combine(cell2mat(idxRegion),:);
                
                % save([folderMaskExtract  filesep 'groupLists.mat'],'groupLists');
            end
            
        end %
    else
        for iExtract = 1:numel(maskParamsCurrent.extractType)
            %             extractType = char(maskParamsCurrent.extractType{iExtract});
            % %
            %             folderMaskExtract = [folderMask filesep extractType];
            
            %
            %             if ~isdir(folderMaskExtract)
            %                 mkdir(folderMaskExtract)
            %             end
            %
            
            % comp different spatial regions between groups
            currentGL = cell(length(btwGrpNames));
            for iWind = 1:maskParamsCurrent.numWindows+1
                for iRegion = 1:3 % h patter
                    for iGroup = 1:length(btwGrpNames)
                        currentGL{iGroup} = groupListsSPreCluster{iExtract}{iGroup}{iRegion}{iWind};% extract
                    end
                    
                    groupLists{iExtract}{iRegion}{iWind}= vertcat(currentGL{:});
                    
                end
                
            end
            
            %save([folderMaskExtract filesep 'groupLists.mat'],'groupLists');
            
        end
    end %
    
    %%
    
    
    
    if append == 1;
        
        iMaskSave = numel(micropatternOutput)+iMask;
    else
        iMaskSave = iMask;
    end
    
    micropatternOutput{iMaskSave}.groupList = groupLists;
    micropatternOutput{iMaskSave}.maskParams = maskParamsCurrent;
    micropatternOutput{iMaskSave}.maskParams.exampleSubRoiSet = exampleSubRoiSet; % save these so can read in an example mask for the colormaps
    micropatternOutput{iMaskSave}.maskParams.exampleRoiYXWholeCell = exampleRoiYXWholeCell;
    micropatternOutput{iMaskSave}.groupListWholeCell = groupList;
    save([saveDir filesep 'micropatternOutput'],'micropatternOutput');
    
    
    
    
    
    
    
end % for each mask condition






end

