function []=TFM_part_2_perfSDC()
% INPUT      area [yTL xTL yBR xBR]
%                           yTL  : y coordinate of the top-left corner
%                           xTL  : x coordinate of the top-left corner
%                           yBR  : y coordinate of the bottom-right corner
%                           xBR  : x coordinate of the bottom-right corner
%            Pass area=[] to manually draw a region of interest

load('fileAndFolderNames.mat')

if ~strcmp(pwd,path_ProjFolder)
    display('Before running this script browse to the FSM project folder')
    return
end

%**************************************************************************
% check the SDC
%**************************************************************************
[flowSDC_FileList]=getFileListFromFolder(path_corrSDC_flow,'flow');
numFlowFilesSDC=length(flowSDC_FileList);

nbins=10;
if numFlowFilesSDC>=3
    k=0;
    for i=[1 round((numFlowFilesSDC+1)/2) numFlowFilesSDC]
        %title(['flow: ',num2str(i)])
        fileStructure=load(flowSDC_FileList{i});
        currentFlowFileSDC=fileStructure.flow;
        subplot(3,3,1+k*3)
        quiver(currentFlowFileSDC(:,2),currentFlowFileSDC(:,1),currentFlowFileSDC(:,4)-currentFlowFileSDC(:,2),currentFlowFileSDC(:,3)-currentFlowFileSDC(:,1))
        set(gca,'Ydir','reverse')
        title('Flow')

        subplot(3,3,2+k*3)
        hist(currentFlowFileSDC(:,4)-currentFlowFileSDC(:,2),nbins)
        title('Histogram of the x-component of the Flow')

        subplot(3,3,3+k*3)
        hist(currentFlowFileSDC(:,3)-currentFlowFileSDC(:,1),nbins)
        title('Histogram of the y-component of the Flow')
        k=k+1;
    end
else
    fileStructure=load(flowSDC_FileList{1});
    currentFlowFileSDC=fileStructure.flow;
    checkFlow(currentFlowFileSDC,nbins)
    title(['flow: ',num2str(1)])
end

reply = input('Is the detected flow OK? Y/N [Y]: ', 's');
if strcmp(reply,'n') || strcmp(reply,'N') || strcmp(reply,'no')
    display('!!!Don t go on with the analysis!!!')
    display('After this script has completed do the following:')
    display(['Copy: ',folder_name_BeadsRegPix,' -> ',folder_name_Beads])
    display(['Copy: ',folder_name_CellsRegSub,' -> ',folder_name_Cells])
    display(['Copy: ',folder_name_XtraRegSub ,' -> ',folder_name_Xtra])
    display(['Copy: ',folder_name_Xtra2ndRegSub ,' -> ',folder_name_Xtra2nd])
    display('Finally delete all created folders')       
    
    reg_method_Cell='pixelwise';
    
    redo_analysis=1;
else
    reg_method_Cell='interpolative';
    redo_analysis=0;
end


% Calculate the StageDriftCorrection:
[T,T_path]=calcStageDriftUsingCorrFlow(flowSDC_FileList,path_mechSDC);


%**************************************************************************
% Make the Bead folders
%**************************************************************************
display('Making the Bead folders')

createRegFolder(path_BeadsCollapsed,path_BeadsRegPix,[],T_path,'pixelwise');

% [sortedBeadsCollapsedFileList]=getFileListFromFolder(path_BeadsCollapsed);
% 
% % Check if the dir exists:
% if ~isdir(path_BeadsRegPix)
%     mkdir(path_BeadsRegPix)
% end
% 
% % Register the Bead Images pixelwise:
% perfRegInPixStep('pixelwise', sortedBeadsCollapsedFileList, path_BeadsRegPix, T_path)

%**************************************************************************
% Make the Cell folders
%**************************************************************************
display('Making the Cell folders')

createRegFolder(path_CellsCollapsed,path_CellsRegSub,path_CellsFinal,T_path,reg_method_Cell);

% [sortedCellsCollapsedFileList]=getFileListFromFolder(path_CellsCollapsed);
% 
% % Check if the dir exists:
% if ~isdir(path_CellsRegSub)
%     mkdir(path_CellsRegSub)
% end
% 
% % Sub-pixel registration of cell images:
% perfRegInPixStep(reg_method_Cell, sortedCellsCollapsedFileList, path_CellsRegSub, T_path)
% 
% 
% [sortedCellsRegSubFileList]=getFileListFromFolder(path_CellsRegSub);
% 
% %remove the reference frames:
% % Check if the dir exists:
% if ~isdir(path_CellsFinal)
%     mkdir(path_CellsFinal)
% end
% 
% removeRefFrame(sortedCellsRegSubFileList,path_CellsFinal)
% rmdir(path_CellsRegSub,'s');


%**************************************************************************
% Make the Xtra folders if necessary
%**************************************************************************
display('Making the Xtra folders')

existXtra=createRegFolder(path_XtraCollapsed,path_XtraRegSub,path_XtraFinal,T_path,reg_method_Cell);

% [sortedXtraCollapsedFileList]=getFileListFromFolder(path_XtraCollapsed);
%
% test_Xtra=~isempty(sortedXtraCollapsedFileList);
% if test_Xtra
%     display('Making the Xtra folders')
%     % Check if the dir exists:
%     if ~isdir(path_XtraRegSub)
%         mkdir(path_XtraRegSub)
%     end
% 
%     % Sub-pixel registration of cell images:
%     perfRegInPixStep(reg_method_Cell, sortedXtraCollapsedFileList, path_XtraRegSub, T_path)
% 
% 
%     [sortedXtraRegSubFileList]=getFileListFromFolder(path_XtraRegSub);
% 
%     %remove the reference frames:
%     % Check if the dir exists:
%     if ~isdir(path_XtraFinal)
%         mkdir(path_XtraFinal)
%     end
% 
%     removeRefFrame(sortedXtraRegSubFileList,path_XtraFinal)
%     % rmdir(path_XtraRegSub,'s'); otherwise error by redoing the analysis
%     % on bad datasets
% end

%**************************************************************************
% Make the Xtra folders if necessary
%**************************************************************************
display('Making the Xtra2nd folders')

existXtra2nd=createRegFolder(path_Xtra2ndCollapsed,path_Xtra2ndRegSub,path_Xtra2ndFinal,T_path,reg_method_Cell);

%**************************************************************************
% Finish
%**************************************************************************

if redo_analysis==1
    reply = input('Do you want the program to perform the file operations (copy and delete folders)? Y/N [N]: ', 's');
    if strcmp(reply,'y') || strcmp(reply,'Y') || strcmp(reply,'yes')
        
        %This should be copied to a backup-data:
        folder_backup=[path_DataFolder,filesep,'backup_original_files'];
        if ~isdir(folder_backup)
            mkdir(folder_backup)
        end        
        
        movefile(path_BeadsFolder   ,folder_backup);
        movefile(path_CellsFolder   ,folder_backup);
        movefile(path_RefFrameFolder,folder_backup);        
        if existXtra
            movefile(path_XtraFolder    ,folder_backup);
        end
        if existXtra2nd
            movefile(path_Xtra2ndFolder    ,folder_backup);
        end
        
        % move one of the padded Ref-frames to the reference Frame Folder
        [sortedBeadsRegSubFileList]=getFileListFromFolder(path_BeadsRegPix);
        if ~isdir(path_RefFrameFolder)
            mkdir(path_RefFrameFolder)
        end
        copyfile(sortedBeadsRegSubFileList{1},path_RefFrameFolder);
        
        % Now remove the reference frames:
        % Check if the dir exists:
        path_BeadsFinal= [path_DataFolder,filesep,'BeadsFinal'];
        if ~isdir(path_BeadsFinal)
            mkdir(path_BeadsFinal)
        end
        removeRefFrame(sortedBeadsRegSubFileList,path_BeadsFinal)
        

        movefile([path_BeadsFinal,filesep,'*'] ,path_BeadsFolder,'f');       
        movefile([path_CellsFinal,filesep,'*'] ,path_CellsFolder,'f');
        if existXtra
            movefile([path_XtraFinal,filesep,'*'] ,path_XtraFolder,'f');  
        end
        if existXtra2nd
            movefile([path_Xtra2ndFinal,filesep,'*'] ,path_Xtra2ndFolder,'f');  
        end   
                
        rmdir(path_BeadsCollapsed,'s');        
        rmdir(path_BeadsRegPix,'s');
        rmdir(path_BeadsFinal,'s');
        rmdir(path_BeadsCropSDC,'s');
       
        rmdir(path_CellsCollapsed,'s');
        %rmdir(path_CellsRegSub,'s');
        rmdir(path_CellsFinal,'s');        
             
        if existXtra
            rmdir(path_XtraCollapsed,'s');
            rmdir(path_XtraRegSub,'s');
            rmdir(path_XtraFinal,'s');
        end
        
        if existXtra2nd
            rmdir(path_Xtra2ndCollapsed,'s');
            rmdir(path_Xtra2ndRegSub,'s');
            rmdir(path_Xtra2ndFinal,'s');
        end
        
        %rmdir(path_corrTFM,'s');
        rmdir(path_corrSDC_flow,'s');
        rmdir(path_mechSDC,'s');
            
        return;
    end
end

% This is now done in TFM_part_1
% %create the geomFile:
% [sortedBeadsRegPixFileList]=getFileListFromFolder(path_BeadsRegPix);
% imageInform=imfinfo(sortedBeadsRegPixFileList{1});
% 
% leftUpperCorner=[1 1];
% lowerRightCorner=[imageInform.Width imageInform.Height];
% edgeErrodeTFM=ceil(max(max(abs(T))))+25;
% 
% if ~isdir(path_corrTFM)
%     mkdir(path_corrTFM)
% end
% 
% [fieldGeom]=createFieldGeom(leftUpperCorner,lowerRightCorner,edgeErrodeTFM,path_corrTFM);
% display('done')