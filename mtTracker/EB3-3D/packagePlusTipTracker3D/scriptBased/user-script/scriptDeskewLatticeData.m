%% load an example of a collection of sequence folders (paths in a matlab Cell) 
goodCells=exampleCellSelection();

%% deskew each sequence described in <goodCells>, save the results in:
% '/project/cellbiology/gdanuser/shared/proudot/project/EB3-3D-track/data-analysis/complete-datasets/goodCells'
%
% When assigned, the <rootPath> '/work/gdanuser/proudot/project/EB3-3D-track/data-analysis/completeDataset'
% is used as a common root for the paths in <goodCells>, i.e. the file tree
% structure is conserved until that point in the file tree. 
%
% If <rootPath> is not provided, the results are saved in a separate
% 'deskew' folder in each original sequence folder. Then, only the movieList file
% is saved in the outputDirectory. 
deskewLatticeData(goodCells,'/project/cellbiology/gdanuser/shared/proudot/project/EB3-3D-track/data-analysis/complete-datasets/goodCellsTest/','rootPath','/home2/proudot/project/EB3-3D-track/data-analysis/completeDataset/')
