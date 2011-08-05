function [forceField]=TFM_part_3_calcForces(method,xrange,yrange,doShift,doRotReg,solMethodBEM)
% INPUT      method: 'FastBEM' or 'FTTC'
%            area [yTL xTL yBR xBR]
%                           yTL  : y coordinate of the top-left corner
%                           xTL  : x coordinate of the top-left corner
%                           yBR  : y coordinate of the bottom-right corner
%                           xBR  : x coordinate of the bottom-right corner
%            Pass area=[] to manually draw a region of interest

if nargin < 1 || isempty(method)
    method ='FastBEM';
    %method='FTTC';
end

if nargin <2 || isempty(xrange) || isempty(yrange)
    xrange=[];
    yrange=[];
end

if nargin <4 || isempty(doShift)
    doShift=true;
end

if nargin <5 || isempty(doRotReg)
    doRotReg=0;
end

if nargin <6 || isempty(solMethodBEM)
    solMethodBEM='QR';
end

load('fileAndFolderNames.mat')

if ~strcmp(pwd,path_ProjFolder)
    display('Before running this script browse to the FSM project folder')
    return
end

%**************************************************************************
% calculate displacement fields
%**************************************************************************
[flowTFM_FileList]=getFileListFromFolder(path_corrTFM_flow,'flow');
numFlowFilesTFM=length(flowTFM_FileList);


display('for default values just press enter:')
display('The Poissons ratio is set to the   default:      0.5   ');

yModu_kPa    = input('Young module in kPa               (default:   20 )   : ');
yModu_Pa = yModu_kPa*1000;

pixSize_mu  = input('Pixel size in mu                  (default:    0.163): ');

if strcmp(method,'FastBEM')
    meshPtsFwdSol = input('Number of mesh pts of fwdSolution (default:     2^12): ');
else
    meshPtsFwdSol =0;
end

fieldsOK='n';
while strcmp(fieldsOK,'n') || strcmp(fieldsOK,'no') || strcmp(fieldsOK,'N')
    filter      = input('Filter spec: [numStd boxSizeLocFac boxSizeGlbFac]=[18 10 6 3] (default: no filter): ');
    regParam    = input('Regularization parameter          (default:  10^(-7)): ');    

    [displField, forceField, ~]=createDisplField(path_ResidualT,flowTFM_FileList,path_mechTFM,filter,yModu_Pa,[],pixSize_mu,regParam,method,meshPtsFwdSol,xrange,yrange,doRotReg,solMethodBEM);
    
    fieldsOK=input('Are you satisfied with the results?: Y/N [Y]: ','s');
end

% get the time point for each frame:
beadImageList=getFileListFromFolder(path_BeadsFolder);
[t, ~, ~, dt_mean, dt_std]=getTimeList(beadImageList);

for i=1:length(displField)
    displField(i).par.t          =t(i);
    displField(i).par.dt_mean    =dt_mean;
    displField(i).par.dt_std     =dt_std;
    
    forceField(i).par.t          =t(i);
    forceField(i).par.dt_mean    =dt_mean;
    forceField(i).par.dt_std     =dt_std;
end    
    

if doShift
    [forceField]=shiftForceField(forceField,displField);
    %save([path_mechTFM,filesep, 'forceField.mat'], 'forceField', '-v7.3');
end

save([path_mechTFM,filesep, 'forceField.mat'], 'forceField', '-v7.3');
save([path_mechTFM,filesep, 'displField.mat'], 'displField', '-v7.3');

return;

sortedCellFinalFileList=getFileListFromFolder(path_CellsFinal);

% Check if the dir exists:
if ~isdir(path_cellsWithDispl)
    mkdir(path_cellsWithDispl)
end

% Check if the dir exists:
if ~isdir(path_cellsWithforce)
    mkdir(path_cellsWithforce)
end


plotCellsWithMech(sortedCellFinalFileList,displField,path_cellsWithDispl,'displField')
plotCellsWithMech(sortedCellFinalFileList,forceField,path_cellsWithforce,'forceField')
