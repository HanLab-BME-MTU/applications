recreatePathAndFolderNames;
load('fileAndFolderNames.mat');

pattern_RefFrame ='w2642';
% get the path to the reference frame:
[refFrameFileList]=   getFileListFromFolder(path_RefFrameFolder, pattern_RefFrame);
% if there is more than one ref. Frame in the List, take the first one:
path_RefFrame=refFrameFileList{1};

createCollapsedFolder(path_Xtra2ndFolder,path_Xtra2ndWithRefFrame,path_Xtra2ndCollapsed,path_RefFrame);

createPreRegFolder(path_Xtra2ndCollapsed,path_tmp,path_PreRegT);

T_path=[path_mechSDC,filesep,'StageDriftCorrectionByCorrFlow.mat'];

reg_method_Cell='interpolative';

existXtra2nd=createRegFolder(path_Xtra2ndCollapsed,path_Xtra2ndRegSub,path_Xtra2ndFinal,T_path,reg_method_Cell);