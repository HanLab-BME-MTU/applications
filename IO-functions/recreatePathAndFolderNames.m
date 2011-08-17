function recreatePathAndFolderNames()

[path_ProjFolder body no ext]=getFilenameBody(pwd);
if ~strcmp(body,'data')
    display('Before running this script browse to the "data" folder')
    return
end
path_DataFolder=[path_ProjFolder,filesep,body];

pattern_WellName='well';
startIndex=strfind(path_ProjFolder, pattern_WellName);
well_name=path_ProjFolder(startIndex:end);

% The folder for these experiments is:
path_experiment=path_ProjFolder(1:startIndex-2); %-2 because there is a filesep



% Set TFM parameters:
pattern_RefFrame                        ='w2642';
folder_name_Beads                       ='Beads';
folder_name_Cells                       ='Cells';
folder_name_Xtra                        ='Ecad'; %This Xtra-Folder (for Caax or Ecad) is treated as the Cell images.
folder_name_Xtra2nd                     ='Nuclei'; %This Xtra-Folder (for Caax or Ecad) is treated as the Cell images.
folder_name_RefFrame                    ='Reference Frame';


folder_name_BeadsWithRefFrame           ='2BeadsWithRefFrame';
folder_name_CellsWithRefFrame           ='2CellsWithRefFrame';
folder_name_XtraWithRefFrame            =['2',folder_name_Xtra,'WithRefFrame'];
folder_name_Xtra2ndWithRefFrame         =['2',folder_name_Xtra2nd,'WithRefFrame'];
folder_name_BeadsCollapsed              ='3BeadsCollapsed';
folder_name_CellsCollapsed              ='3CellsCollapsed';
folder_name_XtraCollapsed               =['3',folder_name_Xtra,'Collapsed'];
folder_name_Xtra2ndCollapsed            =['3',folder_name_Xtra2nd,'Collapsed'];
folder_name_BeadsCropSDC                ='4BeadsCropSDC';
folder_name_tmp                         ='PreRegTmpFolder';

% These folder names are prepared for TFM_part_2_perfSDC.m only, they are not
% used here yet:
folder_name_BeadsRegPix                 ='5BeadsRegPix';
folder_name_CellsRegSub                 ='5CellsRegSub';
folder_name_XtraRegSub                  =['5',folder_name_Xtra,'RegSub'];
folder_name_Xtra2ndRegSub               =['5',folder_name_Xtra2nd,'RegSub'];
folder_name_CellsFinal                  ='6CellsFinal';
folder_name_XtraFinal                   =['6',folder_name_Xtra,'Final'];
folder_name_Xtra2ndFinal                =['6',folder_name_Xtra2nd,'Final'];

folder_name_corrSDC_flow                =['corrSDC',filesep,'flow'];
folder_name_corrTFM                     ='corrTFM';
folder_name_corrTFM_flow                =['corrTFM',filesep,'flow'];
folder_name_mechSDC                     ='mechSDC';
folder_name_mechTFM                     ='mechTFM';
folder_name_displFieldWithCells         ='displFieldWithCells';
folder_name_forceFieldWithCells         ='forceFieldWithCells';
folder_name_result                      ='result';

path_BeadsFolder      =[path_DataFolder,filesep,folder_name_Beads];
path_CellsFolder      =[path_DataFolder,filesep,folder_name_Cells];
path_XtraFolder       =[path_DataFolder,filesep,folder_name_Xtra];
path_Xtra2ndFolder    =[path_DataFolder,filesep,folder_name_Xtra2nd];
path_RefFrameFolder   =[path_DataFolder,filesep,folder_name_RefFrame];
path_BeadsWithRefFrame=[path_DataFolder,filesep,folder_name_BeadsWithRefFrame];
path_CellsWithRefFrame=[path_DataFolder,filesep,folder_name_CellsWithRefFrame];
path_XtraWithRefFrame =[path_DataFolder,filesep,folder_name_XtraWithRefFrame];
path_Xtra2ndWithRefFrame =[path_DataFolder,filesep,folder_name_Xtra2ndWithRefFrame];
path_BeadsCollapsed   =[path_DataFolder,filesep,folder_name_BeadsCollapsed];
path_CellsCollapsed   =[path_DataFolder,filesep,folder_name_CellsCollapsed];
path_XtraCollapsed    =[path_DataFolder,filesep,folder_name_XtraCollapsed];
path_Xtra2ndCollapsed =[path_DataFolder,filesep,folder_name_Xtra2ndCollapsed];
path_BeadsCropSDC     =[path_DataFolder,filesep,folder_name_BeadsCropSDC];
path_tmp              =[path_DataFolder,filesep,folder_name_tmp];

% These pathnames are prepared for TFM_part_2_perfSDC.m only, they are not
% used here yet:
path_BeadsRegPix   =[path_DataFolder,filesep,folder_name_BeadsRegPix];
path_CellsRegSub   =[path_DataFolder,filesep,folder_name_CellsRegSub];
path_XtraRegSub    =[path_DataFolder,filesep,folder_name_XtraRegSub];
path_Xtra2ndRegSub =[path_DataFolder,filesep,folder_name_Xtra2ndRegSub];
path_CellsFinal    =[path_DataFolder,filesep,folder_name_CellsFinal];
path_XtraFinal     =[path_DataFolder,filesep,folder_name_XtraFinal];
path_Xtra2ndFinal  =[path_DataFolder,filesep,folder_name_Xtra2ndFinal];

path_corrSDC_flow  =[path_ProjFolder,filesep,folder_name_corrSDC_flow];
path_corrTFM       =[path_ProjFolder,filesep,folder_name_corrTFM];
path_corrTFM_flow  =[path_ProjFolder,filesep,folder_name_corrTFM_flow];
path_mechSDC       =[path_ProjFolder,filesep,folder_name_mechSDC];
path_mechTFM       =[path_ProjFolder,filesep,folder_name_mechTFM];
path_cellsWithDispl=[path_mechTFM   ,filesep,folder_name_displFieldWithCells];
path_cellsWithforce=[path_mechTFM   ,filesep,folder_name_forceFieldWithCells];
path_result_dir    =[path_experiment,filesep,folder_name_result];

% These are the result file names:
file_name_PreRegT        = 'PreRegByXCorrT.mat';
file_name_ResidualT      = 'ResdualT.mat';
file_name_displField     = 'displField.mat';
file_name_forceField     = 'forceField.mat';
file_name_cellCellForces = 'cellCellForces.mat';
file_name_trackedNet     = 'trackedNet.mat';

% These is the path to the pre-registration folder:
path_PreRegT = [path_DataFolder,filesep,file_name_PreRegT];

% These pathnames are prepared for TFM_part_3_calcForces:
path_ResidualT  =[path_mechSDC,filesep,file_name_ResidualT];


% These pathnames are prepared for TFM_part_4_clusterAnalysis:
path_displField    =[path_mechTFM,filesep,file_name_displField];
path_forceField    =[path_mechTFM,filesep,file_name_forceField];
path_cellCellForces=[path_mechTFM,filesep,file_name_cellCellForces];
path_trackedNet    =[path_mechTFM,filesep,file_name_trackedNet];

cd ..

save('fileAndFolderNames.mat','folder*','path*','file_name*','well_name')