function [pathAnalysisAll, MLNames, groupNames, usedSelectedFoldersMat,...
    specificName,refDirTif, MLdirect]=chooseSelectedFolders
%function [pathAnalysisAll, MLNames, groupNames,
%usedSelectedFoldersMat]=chooseSelectedFolders let users choose predefined
%selectedfolder.mat or directly choose ML files to perform many batch
%functions such as strainEnergyBatch, adhesionAnalysisBatch,
%postAdhAnalysisPackBatch or createSingleFrameMDsFromMultiFrameMDs.
% Sangyoon Han Aug 2020

usedSelectedFoldersMat=false;
isDesktopAvail = usejava('desktop');

if isDesktopAvail
    [fileSFolders, pathSFolders] = uigetfile('*.mat','Select selectedFolders.mat.  If do not have one, click cancel');
else
    disp({'Type selectedFolders.mat.  If do not have one, push enter';
        ['Your current path: ' pwd]});
    rawPath = input(': ','s');
    if isempty(rawPath)
        pathSFolders = 0;
    else
        [pathSFolders, fileSFolders] = fileparts(rawPath);
    end
end

groupNames=[];
MLdirect=false;
refDirTif=[];

if ~ischar(pathSFolders) && pathSFolders==0
    analysisFolderSelectionDone = false;
    ii=0;
    rootFolder=pwd;
    while ~analysisFolderSelectionDone
        ii=ii+1;
        if isDesktopAvail
            curPathProject = uigetdir(rootFolder,'Select each analysis folder that contains movieList.mat (Click Cancel when no more)');
        else
            disp({'Type each analysis folder that contains movieList.mat.  If do not have one, push enter';
                ['Your current path: ' pwd]});
            rawPath = input(': ','s');
            if isempty(rawPath)
                curPathProject = 0;
            else
                curPathProject = rawPath;
            end
        end
        if ~ischar(curPathProject) && curPathProject==0
            analysisFolderSelectionDone=true;
            MLdirectNeeded=true;
        else
%             analysisFolderSelectionDone=true;
            [curPathProject2,finalFolder] = fileparts(curPathProject);
            pathAnalysisAll{ii} = curPathProject;
            if isempty(finalFolder)
                [~,finalFolder] = fileparts(curPathProject2);
            end
            groupNames{ii} = finalFolder;
            MLNames{ii} = 'movieList.mat';
            MLdirectNeeded=false;
        end
    end
    if ii==1 && MLdirectNeeded
        MLSelectionDone = false;
        ii=0;
        while ~MLSelectionDone
            ii=ii+1;
            if isDesktopAvail
                [nameML, curPathML] = uigetfile('*.mat','Select each movieList (Click Cancel when no more)');
            else
                disp({'Select each movieList.  If do not have one, push enter';
                    ['Your current path: ' pwd]});
                rawPath = input(': ','s');
                if isempty(rawPath)
                    curPathML = rawPath;
                else
                    [curPathML, nameML] = fileparts(rawPath);
                end
            end
            if ~ischar(curPathML) || isempty(curPathML)
                MLSelectionDone=true;
            else
                curML = load(fullfile(curPathML,nameML),'ML'); curML=curML.ML;
                pathAnalysisAll{ii} = curML.getPath;
                try
                    groupNames{ii} = nameML(10:end-4); %Excluding first 'movieList' and last '.mat'
                catch % when movieList is just movieList.mat, use the name of the containing folder
                    [~,finalFolder] = fileparts(pathAnalysisAll{ii});
                    groupNames{ii} = finalFolder;
                end
                MLNames{ii} = nameML;
            end
        end
        MLdirect=true;
    end
    if isempty(groupNames)
        disp('Nothing is happening')
        return
    end
    specificName = strjoin(groupNames);
    rootAnalysis = pathAnalysisAll{1};
    refDirTif = [];
    save([rootAnalysis filesep 'selectedFolders' specificName '.mat'], 'rootAnalysis','pathAnalysisAll','MLNames','groupNames')
else
    usedSelectedFoldersMat=true;    
    selectedFolders=load([pathSFolders filesep fileSFolders]);
    pathAnalysisAll=selectedFolders.pathAnalysisAll;
    specificName=fileSFolders(16:end-4);
    groupNames=selectedFolders.groupNames;
    if isfield(selectedFolders,'refDirTifAll')
        refDirTif = selectedFolders.refDirTifAll;
    else
        disp('refDirTifAll was not saved right in setupMovieDataFromND!')
    end
    
    try
        MLNames = selectedFolders.MLNames;%'movieList.mat';
    catch
        for k=1:numel(pathAnalysisAll)
            MLNames{k} = 'movieList.mat';
        end
        
    end
end

