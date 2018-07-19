function [constrForceField]=TFM_part_4_cutOutForceField(ROI)
load('fileAndFolderNames.mat');
fileStruct=load(path_forceField);
forceField=fileStruct.forceField;
close all;

if nargin<1
    ROI=[];
end

%**************************************************************************
% Cut out the force field:
%**************************************************************************
[sortedXtraCollapsedFileList]=getFileListFromFolder(path_XtraCollapsed);
test_Xtra=~isempty(sortedXtraCollapsedFileList);

if test_Xtra
    [constrForceField]=cutOutForceFieldManyCells(forceField,path_XtraFinal,path_cellCellForces,[],ROI);
%     % and save to the result dir:    
%     if ~isdir(path_result_dir)
%         mkdir(path_result_dir)
%     end
%     save([path_result_dir,filesep, 'analysisConstrCellCluster_',well_name,'.mat'],'constrForceField','-v7.3');
    
    % plotCellCellForces(constrForceField,forceField,imageFileList,target_dir,doSave)
else
    [constrForceField, constrDisplField]=cutOutForceFieldSingleCell(displField,forceField,path_CellsFinal ,path_mechTFM);
    % and save to the result dir:    
    if ~isdir(path_result_dir)
        mkdir(path_result_dir)
    end
    save([path_result_dir,filesep, 'analysisConstrSingleCell_',well_name,'.mat'],'constrForceField','constrDisplField');
end