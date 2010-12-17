function [constrForceField]=TFM_part_6_clusterAnalysis(constrForceField,forceField)
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

% load the image file list for Ecad or Cytoplasmic marker:
[sortedXtraFinalFileList]=getFileListFromFolder(path_XtraFinal);



%**************************************************************************
% Now perform the Network Analysis:                                       *
%**************************************************************************

% check which frames have been analyzed so far, and which frames contain
% more than one cell, that is, has at least one interface:
toDoList=[];
for frame=1:length(constrForceField)
    if ~isempty(constrForceField{frame}) && isfield(constrForceField{frame},'interface')
        toDoList=horzcat(toDoList,frame);
    end
end

% Identify cells that are special: Radius [100] about the center of mass of
% the cells. Channel and threshold has to be defined by the user:
constrForceField=identifySpecialCells(constrForceField);
    

% perform the analysis on all these frames:
for frame=toDoList    
    constrForceField=perfNetworkAnalysis(constrForceField,frame);
    figure(6)
    plotCellNetwork(constrForceField{frame}.network);
    title(['Network for frame: ',num2str(frame)])
end


for frame=toDoList
    display(['Working on frame: ',num2str(frame)])
    % Now perform the Cluster Analysis:
    tic;
    constrForceField=perfClusterAnalysis(constrForceField,forceField,frame);
    toc;
    plotClusterAnalysis(constrForceField,forceField,sortedXtraFinalFileList,frame,0);
    
    % Replot the Network results:
    figure(7)
    plotCellNetwork(constrForceField{frame}.network);        
end
save(path_cellCellForces, 'constrForceField','-v7.3');



