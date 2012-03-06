function [output] = makeGroupListsMicropatternWithinGroup(projList,maskParams,extractType)
% output: cell array {number of region types} {window Number} Ideal for
% combining to make all lists 


%make all within groupLists
% note there are probably easier ways of clustering 
% but I kept this because it was my first version and it was finished.
% and yeah i know my sideline is filled with orange :0 f-off

     %filepushing: god i hate it 
     upOne = getFilenameBody(char(projList(1))); 
     upTwo = getFilenameBody(upOne); 
     s1 = num2str(maskParams.numWindows); 
     s2 = num2str(maskParams.windowSize); 
     if maskParams.subRegions == 1
         
     s3 = 'subRegions'; 
     else 
         s3 = ''; 
     end 
     s4  = extractType; 
     subRoiFolderName = ['SUBROIS_' s1 '_' s2 '_umWindows_' s3 '_' s4];    
     groupListDir = [upTwo filesep 'ANALYSIS' s1 '_' s2 '_umWindows_' s3 '_' s4];  
     if ~isdir(groupListDir)
         mkdir(groupListDir)
     end 
     
     
     [~, groupName]= getFilenameBody(upTwo); 
  
      numWindows = maskParams.numWindows;
      windowSize = maskParams.windowSize; 
    %%%%% GroupLists Non-Adhesion %%%%%%%%%%
      
     
    %Individual groupLists for Small Windows From Cell Edge
    for iWindow = 1:numWindows
    
     %NonAdDir = [statDir filesep 'NonAdhesion_' num2str(iWindow*windowSize) 'uM'];
     %mkdir(NonAdDir);
   
    groupList = []; 
    groupListCount = 1;
    for iProj = 1:length(projList)
        
        groupList{groupListCount,1} = [ groupName 'NonAdhesion_' num2str(iWindow*windowSize) 'uM'];
        groupList{groupListCount+1,1} = [ groupName 'NonAdhesion_' num2str(iWindow*windowSize) 'uM'];
        
        groupList{groupListCount,2} = [(char(projList(iProj))) filesep subRoiFolderName filesep 'sub_' num2str(iWindow)];
        groupList{groupListCount+1,2} = [(char(projList(iProj))) filesep subRoiFolderName filesep 'sub_' num2str(iWindow+numWindows+1)];
    
        groupListCount = groupListCount + 2;
    end % if Proj
    
    %Create GroupList to Compare Among Windows in Same Region Type
   
    if iWindow == 1
        groupListNonAdCompare = groupList;
    else 
        groupListNonAdCompare = [groupListNonAdCompare; groupList];
    end
    
    %Save GroupList Files for Each Window Size/Region Type 
    save([groupListDir filesep 'groupListNonAd_' num2str(iWindow*windowSize) 'uM'],'groupList');
    output{1}{iWindow} = groupList; 
    % Hold Window Groups for Between Region Comparisons
    groupListNonAd(:,:,iWindow) = groupList;
    
    end  % end iWindow   
    
    % Create GroupList for Central Region
        groupListCount = 1;
        clear groupList
        for iProj= 1:length(projList)
            groupList{groupListCount,1} = [ groupName 'NonAdhesion_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            groupList{groupListCount+1,1} = [ groupName 'NonAdhesion_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            groupList{groupListCount,2} = [char(projList(iProj)) filesep subRoiFolderName filesep 'sub_' num2str(numWindows+1)];
            groupList{groupListCount+1,2} = [char(projList(iProj)) filesep subRoiFolderName filesep 'sub_' num2str(2*(numWindows+1))];
            groupListCount = groupListCount + 2;
        end
        
     %NonAdDir = [statDir filesep 'NonAdhesion_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
     %mkdir(NonAdDir);
     
     %Save Individual GroupList For Central Region
     save([groupListDir filesep 'groupListNonAd_GreaterThan_' num2str(numWindows*windowSize) 'uM'],'groupList');   
     output{1}{numWindows+1} = groupList; 
  
     % Hold Central Region Group for Between Region Comparisons
     groupListNonAd(:,:,numWindows+1) = groupList;
   
     %Add Central Region to Window Comparison GroupList 
     groupListNonAdCompare = [groupListNonAdCompare; groupList];
     
     %Save Window Comparison GroupList 
     groupList = groupListNonAdCompare; 
     save([groupListDir filesep 'groupListCompareWindows_NonAd'],'groupList');
     
     %%%%%% GroupLists Adhesion %%%%%% 
     
       
        clear groupList
       
        %Individual groupLists for Small Windows From Cell Edge
        for iWindow = 1:numWindows
    
        %AdDir = [statDir filesep 'Adhesion_' num2str(iWindow*windowSize) 'uM'];
        %mkdir(AdDir);
   
    
        groupListCount = 1;
        for iProj = 1:length(projList)
            groupList{groupListCount,1} = [ groupName 'Adhesion_' num2str(iWindow*windowSize) 'uM']; 
            groupList{groupListCount+1,1} = [ groupName 'Adhesion_' num2str(iWindow*windowSize) 'uM'];
            windowNumRight = 2*(numWindows+1) + iWindow;
            windowNumLeft = windowNumRight+numWindows + 1;
            groupList{groupListCount,2} = [char(projList(iProj)) filesep subRoiFolderName filesep 'sub_' num2str(windowNumRight)];
            groupList{groupListCount+1,2} = [char(projList(iProj)) filesep subRoiFolderName filesep 'sub_' num2str(windowNumLeft)];
             
            groupListCount = groupListCount +2;
        end  
        
         %Create GroupList to Compare Among Windows in Same Region Type
        if iWindow == 1
            groupListAdCompare = groupList;
        else 
            groupListAdCompare = [groupListAdCompare; groupList];
        end
        
         % Hold Window Groups for Between Region Comparisons
            groupListAd(:,:,iWindow) = groupList;
        
        %Save GroupList Files for Each Window Size/Region Type 
        save([groupListDir filesep 'groupListAd_' num2str(iWindow*windowSize) 'uM'],'groupList');
        %Ad = 1
        output{2}{iWindow} = groupList; 
       
        
       
        
        end %for iWindows
        
       
        
        % Create GroupList for Central Region
        groupListCount = 1;
        clear groupList
        for iProj= 1:length(projList)
            
            groupList{groupListCount,1} = [groupName 'Adhesion_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            groupList{groupListCount+1,1} = [ groupName 'Adhesion_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            windowNumRight = 2*(numWindows+1) + numWindows+1;
            windowNumLeft = windowNumRight+numWindows+1;
            groupList{groupListCount,2} = [char(projList(iProj)) filesep subRoiFolderName filesep 'sub_' num2str(windowNumRight)]; % 
            groupList{groupListCount+1,2} = [char(projList(iProj)) filesep subRoiFolderName filesep 'sub_' num2str(windowNumLeft)];
            groupListCount = groupListCount + 2;
        end
        
     %AdDir = [statDir filesep 'Adhesion_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
     %mkdir(AdDir);
     
     
     %Save Individual GroupList For Central Region
     save([groupListDir filesep 'groupListAd_GreaterThan_' num2str(numWindows*windowSize) 'uM'],'groupList');   
     output{2}{numWindows+1} = groupList; 
     
     % Hold Central Region Group for Between Region Comparisons
     groupListAd(:,:,numWindows+1) = groupList;
     
     %Add Central Region to Window Comparison GroupList 
     groupListAdCompare = [groupListAdCompare; groupList];
     
     %Save Window Comparison GroupList 
     groupList = groupListAdCompare;
     save([groupListDir filesep 'groupListCompareWindows_Ad'],'groupList');
     
 %%    
     
      %%%%%% GroupLists Adhesion Corners %%%%%% 
      clear groupList
       

    for iWindow = 1:numWindows
    
    % AdCornDir = [statDir filesep 'AdhesionCorn_' num2str(iWindow*windowSize) 'uM'];
    % mkdir(AdCornDir);
   
   
    groupListCount = 1;
    
    for iProj = 1:length(projList)
        
        groupList{groupListCount,1} = [ groupName 'AdhesionCorn_' num2str(iWindow*windowSize) 'uM'];
        groupList{groupListCount+1,1} = [ groupName 'AdhesionCorn_' num2str(iWindow*windowSize) 'uM'];
        groupList{groupListCount+2,1} = [ groupName 'AdhesionCorn_' num2str(iWindow*windowSize) 'uM'];
        groupList{groupListCount+3,1} = [ groupName 'AdhesionCorn_' num2str(iWindow*windowSize) 'uM'];
        
        windowNum1 = 4*(numWindows+1) + iWindow;
        windowNum2 = windowNum1 + numWindows +1;
        windowNum3 = windowNum2 +numWindows +1;
        windowNum4 = windowNum3 + numWindows +1;
        
        groupList{groupListCount,2} = [char(projList(iProj)) filesep subRoiFolderName filesep 'sub_' num2str(windowNum1)];
        groupList{groupListCount+1,2} = [char(projList(iProj)) filesep subRoiFolderName filesep 'sub_' num2str(windowNum2)];
        groupList{groupListCount+2,2} = [char(projList(iProj)) filesep subRoiFolderName filesep 'sub_' num2str(windowNum3)];
        groupList{groupListCount+3,2} = [char(projList(iProj)) filesep subRoiFolderName filesep 'sub_' num2str(windowNum4)];
    
        groupListCount = groupListCount + 4;
           
    end % for iProj
    
    
    %Create GroupList to Compare Among Windows in Same Region Type
    if iWindow == 1
        groupListAdCornCompare = groupList;
    else 
        groupListAdCornCompare = [groupListAdCornCompare; groupList];
    end % end if
     
    % Hold Window Groups for Between Region Comparisons
            groupListAdCorn(:,:,iWindow) = groupList;
    output{3}{iWindow} = groupList; 
    %Save GroupList Files for Each Window Size/Region Type 
    save([groupListDir filesep 'groupListAdCorn_' num2str(iWindow*windowSize) 'uM'],'groupList');
    
    end   % for iWindows 
 
    
     %Create GroupList for Central Region
        groupListCount = 1;
        clear groupList
        for iProj= 1:length(projList)
            
            groupList{groupListCount,1} = ['grp_' groupName 'AdhesionCorn_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            groupList{groupListCount+1,1} = ['grp_' groupName 'AdhesionCorn_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            groupList{groupListCount+2,1} = ['grp_' groupName 'AdhesionCorn_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            groupList{groupListCount+3,1} = ['grp_' groupName 'AdhesionCorn_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            
            windowNum1 = 4*(numWindows+1) + numWindows +1;
            windowNum2 = windowNum1 + numWindows + 1;
            windowNum3 = windowNum2 + numWindows +1;
            windowNum4 = windowNum3 + numWindows +1;
            
            groupList{groupListCount,2} = [char(projList(iProj)) filesep subRoiFolderName filesep 'sub_' num2str(windowNum1)]; 
            groupList{groupListCount+1,2} = [char(projList(iProj)) filesep subRoiFolderName filesep 'sub_' num2str(windowNum2)];
            groupList{groupListCount+2,2} = [char(projList(iProj)) filesep subRoiFolderName filesep 'sub_' num2str(windowNum3)];
            groupList{groupListCount+3,2} = [char(projList(iProj)) filesep subRoiFolderName filesep 'sub_' num2str(windowNum4)];
            
            groupListCount = groupListCount + 4;
        end % for iProj
        
    % AdCornDir = [statDir filesep 'AdhesionCorn_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
    % mkdir(AdCornDir);
    
    %Save Individual Group List for Central Region 
    save([groupListDir filesep 'groupListAdCorn_GreaterThan_' num2str(numWindows*windowSize) 'uM'],'groupList');   
    output{3}{numWindows+1}= groupList; 
    
    % Hold Central Region Group for Between Region Comparisons
    groupListAdCorn(:,:,numWindows+1) = groupList;
    
    %Add Central Region to Window Comparison GroupList 
     groupListAdCornCompare = [groupListAdCornCompare; groupList];
     
     %Save Window Comparison GroupList 
     groupList = groupListAdCornCompare;
     save([groupListDir filesep 'groupListCompareWindows_AdCorn'],'groupList');
     
     
    
    
     clear groupList
     
 %%%%%%%%%%%%%%%%%% GroupList Same Window Different Region Type %%%%%%%
    for iWindow = 1:numWindows
        groupList = [groupListAd(:,:,iWindow) ; groupListAdCorn(:,:,iWindow); groupListNonAd(:,:,iWindow)];
        save([groupListDir filesep 'groupListCompareBtwRegionTypes_' num2str(iWindow*windowSize) 'uM'],'groupList');
    end % for iWindow
    
    groupList = [groupListAd(:,:,numWindows+1) ; groupListAdCorn(:,:,numWindows+1) ; groupListNonAd(:,:,numWindows+1)];
    save([groupListDir filesep 'groupListCompareBtwRegionTypes_GreaterThan' num2str(numWindows*windowSize) 'uM'],'groupList');
    
   %[groupData]  = plusTipPoolGroupData(groupList,1,0,1,1);
   
  %[discrimMat] =  plusTipTestDistr(groupData,[]); 
   %[ ] = ('stats',groupData,[],1,1,1,0);
     
     
    %%%%% GroupLists Non-Adhesion %%%%%%%%%%
    
    %Individual groupLists for Small Windows From Cell Edge
    for iWindow = 1:numWindows
    
     %NonAdDir = [statDir filesep 'NonAdhesion_' num2str(iWindow*windowSize) 'uM'];
     %mkdir(NonAdDir);
    
    groupList = []; 
    groupListCount = 1;
    for iProj = 1:length(projList)
        
        groupList{groupListCount,1} = [ groupName 'NonAdhesion_' num2str(iWindow*windowSize) 'uM'];
        groupList{groupListCount+1,1} = [ groupName 'NonAdhesion_' num2str(iWindow*windowSize) 'uM'];
        
        groupList{groupListCount,2} = [char(projList{iProj}) filesep subRoiFolderName filesep 'sub_' num2str(iWindow)];
        groupList{groupListCount+1,2} = [char(projList{iProj}) filesep subRoiFolderName filesep 'sub_' num2str(iWindow+numWindows+1)];
    
        groupListCount = groupListCount + 2;
    end % if Proj
    
    
    %Create GroupList to Compare Among Windows in Same Region Type
    if iWindow == 1
        groupListNonAdCompare = groupList;
    else 
        groupListNonAdCompare = [groupListNonAdCompare; groupList];
    end
    
    %Save GroupList Files for Each Window Size/Region Type 
    save([groupListDir filesep 'groupListNonAd_' num2str(iWindow*windowSize) 'uM'],'groupList');
    output{1}{iWindow} = groupList; 
    
    % Hold Window Groups for Between Region Comparisons
    groupListNonAd(:,:,iWindow) = groupList;
    
    end  % end iWindow   
    
    % Create GroupList for Central Region
        groupListCount = 1;
        clear groupList
        for iProj= 1:length(projList)
            groupList{groupListCount,1} = [ groupName 'NonAdhesion_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            groupList{groupListCount+1,1} = [ groupName 'NonAdhesion_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            groupList{groupListCount,2} = [char(projList{iProj}) filesep subRoiFolderName filesep 'sub_' num2str(numWindows+1)];
            groupList{groupListCount+1,2} = [char(projList{iProj}) filesep subRoiFolderName filesep 'sub_' num2str(2*(numWindows+1))];
            groupListCount = groupListCount + 2;
        end
        
     %NonAdDir = [statDir filesep 'NonAdhesion_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
     %mkdir(NonAdDir);
     
     %Save Individual GroupList For Central Region
     save([groupListDir filesep 'groupListNonAd_GreaterThan_' num2str(numWindows*windowSize) 'uM'],'groupList');   
     
     % Hold Central Region Group for Between Region Comparisons
     groupListNonAd(:,:,numWindows+1) = groupList;
   
     %Add Central Region to Window Comparison GroupList 
     groupListNonAdCompare = [groupListNonAdCompare; groupList];
     
     %Save Window Comparison GroupList 
     groupList = groupListNonAdCompare; 
     save([groupListDir filesep 'groupListCompareWindows_NonAd'],'groupList');
     
     %%%%%% GroupLists Adhesion %%%%%% 
     
       
        clear groupList
       groupList = [];
        %Individual groupLists for Small Windows From Cell Edge
        for iWindow = 1:numWindows
    
        %AdDir = [statDir filesep 'Adhesion_' num2str(iWindow*windowSize) 'uM'];
        %mkdir(AdDir);
   
    
        groupListCount = 1;
        for iProj = 1:length(projList)
            groupList{groupListCount,1} = [ groupName 'Adhesion_' num2str(iWindow*windowSize) 'uM']; 
            groupList{groupListCount+1,1} = [ groupName 'Adhesion_' num2str(iWindow*windowSize) 'uM'];
            windowNumRight = 2*(numWindows+1) + iWindow;
            windowNumLeft = windowNumRight+numWindows + 1;
            groupList{groupListCount,2} = [char(projList{iProj}) filesep subRoiFolderName filesep 'sub_' num2str(windowNumRight)];
            groupList{groupListCount+1,2} = [char(projList{iProj}) filesep subRoiFolderName filesep 'sub_' num2str(windowNumLeft)];
             
            groupListCount = groupListCount +2;
        end  
        
         %Create GroupList to Compare Among Windows in Same Region Type
        if iWindow == 1
            groupListAdCompare = groupList;
        else 
            groupListAdCompare = [groupListAdCompare; groupList];
        end
        
         % Hold Window Groups for Between Region Comparisons
            groupListAd(:,:,iWindow) = groupList;
        
        %Save GroupList Files for Each Window Size/Region Type 
        save([groupListDir filesep 'groupListAd_' num2str(iWindow*windowSize) 'uM'],'groupList');
        
        
       
        
        end %for iWindows
        
       
        
        % Create GroupList for Central Region
        groupListCount = 1;
        clear groupList
        groupList = []; 
        for iProj= 1:length(projList)
            
            groupList{groupListCount,1} = [groupName 'Adhesion_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            groupList{groupListCount+1,1} = [ groupName 'Adhesion_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            windowNumRight = 2*(numWindows+1) + numWindows+1;
            windowNumLeft = windowNumRight+numWindows+1;
            groupList{groupListCount,2} = [char(projList{iProj}) filesep subRoiFolderName filesep 'sub_' num2str(windowNumRight)]; % 
            groupList{groupListCount+1,2} = [char(projList{iProj}) filesep subRoiFolderName filesep 'sub_' num2str(windowNumLeft)];
            groupListCount = groupListCount + 2;
        end
        
     %AdDir = [statDir filesep 'Adhesion_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
     %mkdir(AdDir);
     
     
     %Save Individual GroupList For Central Region
     save([groupListDir filesep 'groupListAd_GreaterThan_' num2str(numWindows*windowSize) 'uM'],'groupList');   
     
     % Hold Central Region Group for Between Region Comparisons
     groupListAd(:,:,numWindows+1) = groupList;
     
     %Add Central Region to Window Comparison GroupList 
     groupListAdCompare = [groupListAdCompare; groupList];
     
     %Save Window Comparison GroupList 
     groupList = groupListAdCompare;
     save([groupListDir filesep 'groupListCompareWindows_Ad'],'groupList');
     
 %%    
     
      %%%%%% GroupLists Adhesion Corners %%%%%% 
      clear groupList
       
groupList = [];
    for iWindow = 1:numWindows
    
    % AdCornDir = [statDir filesep 'AdhesionCorn_' num2str(iWindow*windowSize) 'uM'];
    % mkdir(AdCornDir);
   
   
    groupListCount = 1;
    
    for iProj = 1:length(projList)
        
        groupList{groupListCount,1} = [ groupName 'AdhesionCorn_' num2str(iWindow*windowSize) 'uM'];
        groupList{groupListCount+1,1} = [ groupName 'AdhesionCorn_' num2str(iWindow*windowSize) 'uM'];
        groupList{groupListCount+2,1} = [ groupName 'AdhesionCorn_' num2str(iWindow*windowSize) 'uM'];
        groupList{groupListCount+3,1} = [ groupName 'AdhesionCorn_' num2str(iWindow*windowSize) 'uM'];
        
        windowNum1 = 4*(numWindows+1) + iWindow;
        windowNum2 = windowNum1 + numWindows +1;
        windowNum3 = windowNum2 +numWindows +1;
        windowNum4 = windowNum3 + numWindows +1;
        
        groupList{groupListCount,2} = [char(projList{iProj}) filesep subRoiFolderName filesep 'sub_' num2str(windowNum1)];
        groupList{groupListCount+1,2} = [char(projList{iProj}) filesep subRoiFolderName filesep 'sub_' num2str(windowNum2)];
        groupList{groupListCount+2,2} = [char(projList{iProj}) filesep subRoiFolderName filesep 'sub_' num2str(windowNum3)];
        groupList{groupListCount+3,2} = [char(projList{iProj}) filesep subRoiFolderName filesep 'sub_' num2str(windowNum4)];
    
        groupListCount = groupListCount + 4;
           
    end % for iProj
    
    
    %Create GroupList to Compare Among Windows in Same Region Type
    if iWindow == 1
        groupListAdCornCompare = groupList;
    else 
        groupListAdCornCompare = [groupListAdCornCompare; groupList];
    end % end if
     
    % Hold Window Groups for Between Region Comparisons
            groupListAdCorn(:,:,iWindow) = groupList;
    
    %Save GroupList Files for Each Window Size/Region Type 
    save([groupListDir filesep 'groupListAdCorn_' num2str(iWindow*windowSize) 'uM'],'groupList');
    
    end   % for iWindows 
 
    
     %Create GroupList for Central Region
        groupListCount = 1;
        clear groupList
        for iProj= 1:length(projList)
            
            groupList{groupListCount,1} = ['grp_' groupName 'AdhesionCorn_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            groupList{groupListCount+1,1} = ['grp_' groupName 'AdhesionCorn_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            groupList{groupListCount+2,1} = ['grp_' groupName 'AdhesionCorn_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            groupList{groupListCount+3,1} = ['grp_' groupName 'AdhesionCorn_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
            
            windowNum1 = 4*(numWindows+1) + numWindows +1;
            windowNum2 = windowNum1 + numWindows + 1;
            windowNum3 = windowNum2 + numWindows +1;
            windowNum4 = windowNum3 + numWindows +1;
            
            groupList{groupListCount,2} = [char(projList{iProj}) filesep subRoiFolderName filesep 'sub_' num2str(windowNum1)]; 
            groupList{groupListCount+1,2} = [char(projList{iProj}) filesep subRoiFolderName filesep 'sub_' num2str(windowNum2)];
            groupList{groupListCount+2,2} = [char(projList{iProj}) filesep subRoiFolderName filesep 'sub_' num2str(windowNum3)];
            groupList{groupListCount+3,2} = [char(projList{iProj}) filesep subRoiFolderName filesep 'sub_' num2str(windowNum4)];
            
            groupListCount = groupListCount + 4;
        end % for iProj
        
    % AdCornDir = [statDir filesep 'AdhesionCorn_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
    % mkdir(AdCornDir);
    
    %Save Individual Group List for Central Region 
    save([groupListDir filesep 'groupListAdCorn_GreaterThan_' num2str(numWindows*windowSize) 'uM'],'groupList');   
    
    % Hold Central Region Group for Between Region Comparisons
    groupListAdCorn(:,:,numWindows+1) = groupList;
    
    %Add Central Region to Window Comparison GroupList 
     groupListAdCornCompare = [groupListAdCornCompare; groupList];
     
     %Save Window Comparison GroupList 
     groupList = groupListAdCornCompare;
     save([groupListDir filesep 'groupListCompareWindows_AdCorn'],'groupList');
     
     
    
    
     clear groupList
     
 %%%%%%%%%%%%%%%%%% GroupList Same Window Different Region Type %%%%%%%
    for iWindow = 1:numWindows
        groupList = [groupListAd(:,:,iWindow) ; groupListAdCorn(:,:,iWindow); groupListNonAd(:,:,iWindow)];
        save([groupListDir filesep 'groupListCompareBtwRegionTypes_' num2str(iWindow*windowSize) 'uM'],'groupList');
    end % for iWindow
    
    groupList = [groupListAd(:,:,numWindows+1) ; groupListAdCorn(:,:,numWindows+1) ; groupListNonAd(:,:,numWindows+1)];
    save([groupListDir filesep 'groupListCompareBtwRegionTypes_GreaterThan' num2str(numWindows*windowSize) 'uM'],'groupList');
    
   %[groupData]  = plusTipPoolGroupData(groupList,1,0,1,1);
   
  %[discrimMat] =  plusTipTestDistr(groupData,[]); 
   %[ ] = ('stats',groupData,[],1,1,1,0);
end