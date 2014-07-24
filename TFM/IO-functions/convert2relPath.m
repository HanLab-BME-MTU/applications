%convert2relPath(path2matFile)
load('fileAndFolderNames.mat');
presentWD=pwd;

if strcmp(pwd,path_ProjFolder)
    display('nothing done')
    return;
end

path_ProjFolder=pwd;

pathNames{1}=path_DataFolder;
pathNames{2}=path_BeadsFolder;
pathNames{3}=path_CellsFolder;
pathNames{4}=path_RefFrameFolder;
pathNames{5}=path_BeadsWithRefFrame;
pathNames{6}=path_CellsWithRefFrame;
pathNames{7}=path_BeadsCollapsed;
pathNames{8}=path_CellsCollapsed;
pathNames{9}=path_BeadsCropSDC;
pathNames{10}=path_BeadsRegPix;
pathNames{11}=path_CellsRegSub;
pathNames{12}=path_CellsFinal;
pathNames{13}=path_corrSDC_flow;
pathNames{14}=path_corrTFM;
pathNames{15}=path_corrTFM_flow;
pathNames{16}=path_mechSDC;
pathNames{17}=path_mechTFM;
pathNames{18}=path_cellsWithDispl;
pathNames{19}=path_cellsWithforce;


% Length of the pattern to be recognized:
pLength=8;
pattern=presentWD(end-pLength:end);

% Find the position of the pattern in the old path
for i=1:length(pathNames)
    posBegin=strfind(pathNames{i},pattern);
    posEnd=posBegin+pLength+2; % 2 since there is also a filesep taking one space

    pathNames{i}=pathNames{i}(posEnd:end);
end



path_DataFolder         =pathNames{1};
path_BeadsFolder        =pathNames{2};
path_CellsFolder        =pathNames{3};
path_RefFrameFolder     =pathNames{4};
path_BeadsWithRefFrame  =pathNames{5};
path_CellsWithRefFrame  =pathNames{6};
path_BeadsCollapsed     =pathNames{7};
path_CellsCollapsed     =pathNames{8};
path_BeadsCropSDC       =pathNames{9};
path_BeadsRegPix        =pathNames{10};
path_CellsRegSub        =pathNames{11};
path_CellsFinal         =pathNames{12};
path_corrSDC_flow       =pathNames{13};
path_corrTFM            =pathNames{14};
path_corrTFM_flow       =pathNames{15};
path_mechSDC            =pathNames{16};
path_mechTFM            =pathNames{17};
path_cellsWithDispl     =pathNames{18};
path_cellsWithforce     =pathNames{19};

clear i pLength pathNames posEnd posBegin pattern presentWD

save('fileAndFolderNames')
