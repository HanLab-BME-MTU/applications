function [constrForceField]=TFM_part_6_clusterAnalysis(constrForceField,forceField)
load('fileAndFolderNames.mat')
doPlot=0;

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



% check which frames have been analyzed so far, and which frames contain
% more than one cell, that is, has at least one interface:
toDoList=[];
for frame=1:length(constrForceField)
    if ~isempty(constrForceField{frame}) && isfield(constrForceField{frame},'interface')
        toDoList=horzcat(toDoList,frame);
    end
end

%**************************************************************************
% 1) Identify myosin cells:                                               *
%**************************************************************************
% Identify cells that are special: Radius [100] about the center of mass of
% the cells. Channel and threshold has to be defined by the user:
imageFileList=getFileListFromFolder(path_Xtra2ndFinal);
constrForceField=identifySpecialCells(constrForceField,imageFileList);
    

%**************************************************************************
% 1) Perform the Network Analysis:                                        *
%**************************************************************************
% perform the analysis on all these frames:
for frame=toDoList    
    constrForceField=perfNetworkAnalysis(constrForceField,frame);
    if doPlot
        figure(6)
        plotCellNetwork(constrForceField{frame}.network);
        title(['Network for frame: ',num2str(frame)])
    end
end


%**************************************************************************
% 2) Measure the intensity along the edges:                               *
%**************************************************************************
imageFileList=getFileListFromFolder(path_XtraFinal);
constrForceField=perfIntMeasures(constrForceField,imageFileList);

%**************************************************************************
% 3) Perform the Cluster Analysis:                                        *
%**************************************************************************
for frame=toDoList
    display(['Working on frame: ',num2str(frame)])
    % Now perform the Cluster Analysis:
    tic;
    constrForceField=perfClusterAnalysis(constrForceField,forceField,frame);
    toc;
    if doPlot
        plotClusterAnalysis(constrForceField,forceField,sortedXtraFinalFileList,frame,0);

        % Replot the Network results:
        figure(7)
        plotCellNetwork(constrForceField{frame}.network);
    end
end

%**************************************************************************
% 4) Track cells and interfaces in the network and create network_tracked *
%**************************************************************************
[constrForceField]=trackNetwork(constrForceField);


%**************************************************************************
% 5) Save the whole data so we can iterate on it later on                 *
%**************************************************************************
save(path_cellCellForces, 'constrForceField','-v7.3');

%**************************************************************************
% 6) Save only the tracked network for statistics (smaller size)          *
%**************************************************************************
for frame=toDoList
    trackedNet{frame}.node =constrForceField{frame}.network_tracked.node;
    trackedNet{frame}.edge =constrForceField{frame}.network_tracked.edge;
    trackedNet{frame}.stats=constrForceField{frame}.network_tracked.stats;
    trackedNet{frame}.par  =constrForceField{frame}.par;
end

[~,fnameFirstBeadImg]=getFileListFromFolder(path_BeadsFolder);
save(path_trackedNet, 'trackedNet','fnameFirstBeadImg','-v7.3');


