function [constrForceField]=TFM_part_5_calcElEnergies(constrForceField,forceField,displField,meshPtsFwdSol)
load('fileAndFolderNames.mat')

if ~strcmp(pwd,path_ProjFolder)
    display('Before running this script browse to the FSM project folder')
    return
end

if nargin < 1 || isempty(constrForceField)
    display('Loading cellCellForces.mat. This file is large and thus takes some while:...')
    tic;
    filestruct=load(path_cellCellForces);
    constrForceField=filestruct.constrForceField;
    toc;    
end

if nargin < 2 || isempty(forceField)
    filestruct=load(path_forceField);
    forceField=filestruct.forceField;
end

if nargin < 3 || isempty(displField)
    filestruct=load(path_displField);
    displField=filestruct.displField;
end

if nargin < 4
    meshPtsFwdSol=input('Enter the number of mesh pts of fwdSolution [2^11]: ');
    if isempty(meshPtsFwdSol)
        meshPtsFwdSol=2^11;
    end
end

% first determine which fields are not empty, then calculate the elastic
% energy for each cell!

for frame=1:length(constrForceField)
    display(['Next frame is: ',num2str(frame)]);
    [constrForceField]=calcElEnergies(constrForceField,forceField,frame,displField,meshPtsFwdSol);
    % here, saving should be still very fast:
end
% or save it after every frame (but that might take some time!):
save(path_cellCellForces, 'constrForceField','-v7.3');