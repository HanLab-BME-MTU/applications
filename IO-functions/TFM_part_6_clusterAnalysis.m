function [constrForceField]=TFM_part_6_clusterAnalysis(constrForceField,forceField,opt,saveAll,doCorrection)
doPlot=0;

if nargin<5
    doCorrection=0;
end

if ~doCorrection
    load('fileAndFolderNames.mat')
end

if ~doCorrection && ~strcmp(pwd,path_ProjFolder)
    display('Before running this script browse to the FSM project folder')
    return;
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

if nargin<3 || isempty(opt)
    opt='all';
end

if nargin<4 || isempty(saveAll)
    saveAll=1;
end

if doPlot && ~doCorrection
    % load the image file list for Ecad or Cytoplasmic marker:
    [sortedXtraFinalFileList]=getFileListFromFolder(path_XtraFinal);
end

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
if ~doCorrection
    imageFileList=getFileListFromFolder(path_Xtra2ndFinal);
    constrForceField=identifySpecialCells(constrForceField,imageFileList,opt);
end
    

%**************************************************************************
% 1) Perform the Network Analysis:                                        *
%**************************************************************************
% toDoList=setdiff(toDoList,[28 29 30]);
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
if doCorrection
    path_XtraFinal='../data/6EcadFinal'; % for the first run I had 'Ecad' the unregistered images.
end
imageFileList=getFileListFromFolder(path_XtraFinal);
constrForceField=perfIntMeasures(constrForceField,imageFileList);

%**************************************************************************
% 3) Perform the Cluster Analysis:                                        *
%**************************************************************************
for frame=toDoList
    display(['Working on frame: ',num2str(frame)])
    % Now perform the Cluster Analysis:
    tic;
    constrForceField=perfClusterAnalysis(constrForceField,forceField,frame,saveAll);
    toc;
    % this plot only works if everything is saved!
    if doPlot && saveAll
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
if ~doCorrection
    save(path_cellCellForces, 'constrForceField','-v7.3');
else
    constrForceFieldCorrected=constrForceField;
    save('cellCellForcesCorrected.mat','constrForceFieldCorrected','-v7.3');
end

%**************************************************************************
% 6) Save only the tracked network for statistics (smaller size)          *
%**************************************************************************
for frame=toDoList
    trackedNet{frame}.node =constrForceField{frame}.network_tracked.node;
    trackedNet{frame}.edge =constrForceField{frame}.network_tracked.edge;
    trackedNet{frame}.stats=constrForceField{frame}.network_tracked.stats;
    trackedNet{frame}.par  =constrForceField{frame}.par;
end

if doCorrection
    path_BeadsFolder='../data/Beads';
end

[~,fnameFirstBeadImg]=getFileListFromFolder(path_BeadsFolder);

if ~doCorrection
    save(path_trackedNet, 'trackedNet','fnameFirstBeadImg','-v7.3');
else
    trackedNetCorrected=trackedNet;
    save('trackedNetCorrected.mat','trackedNetCorrected','fnameFirstBeadImg','-v7.3');
end