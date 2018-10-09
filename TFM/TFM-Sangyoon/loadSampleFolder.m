function [idGroupSelected,MD,curTracksNA]=loadSampleFolder(sampleFolder)
%function [idGroupSelected,MD]=loadSampleFolder(sampleFolder) loads
%idGroupSelected and movieData from the path recorded in sampleFolder.
% Sangyoon Han Sep 2018
    idGroups = load(sampleFolder);
    PathName=fileparts(sampleFolder);
    idGroup1Selected = idGroups.idGroup1Selected;
    idGroup2Selected = idGroups.idGroup2Selected;
    idGroup3Selected = idGroups.idGroup3Selected;
    idGroup4Selected = idGroups.idGroup4Selected;
    idGroup5Selected = idGroups.idGroup5Selected;
    idGroup6Selected = idGroups.idGroup6Selected;
    idGroup7Selected = idGroups.idGroup7Selected;
    idGroup8Selected = idGroups.idGroup8Selected;
    idGroup9Selected = idGroups.idGroup9Selected;
    idGroupSelected={idGroup1Selected,idGroup2Selected,idGroup3Selected,idGroup4Selected,idGroup5Selected,idGroup6Selected,....
                                idGroup7Selected,idGroup8Selected,idGroup9Selected};
    idGroupSelectedTogether = [idGroup1Selected,idGroup2Selected,idGroup3Selected,idGroup4Selected,idGroup5Selected,idGroup6Selected,....
                                idGroup7Selected,idGroup8Selected,idGroup9Selected];
    ss = 0;
    idGroupSelectedReordered = cell(1,9);
    for kk=1:9
        idGroupSelectedReordered{kk} = ss+1:ss+numel(idGroupSelected{kk});
        ss = ss + numel(idGroupSelected{kk});
    end
    try
        mdPath = fileparts(fileparts(PathName));
        curMD = MovieData.load([mdPath filesep 'movieData.mat']);
        curFAPack = curMD.getPackage(curMD.getPackageIndex('FocalAdhesionPackage'));
        curAdhAnalProc = curFAPack.processes_{7};
        iChan=find(curAdhAnalProc.checkChannelOutput,1);
        curTracksNA = curAdhAnalProc.loadChannelOutput(iChan,1,'output','tracksNA','idSelected',idGroupSelectedTogether);
        MD = curMD;
        idGroupSelected = idGroupSelectedReordered;
    catch % backward compatibility
        try
            curImportFilePathTracks = fullfile(PathName,'idsClassified_org.mat');
            curTracksNAfile = load(curImportFilePathTracks,'tracksNA');
            curTracksNA = curTracksNAfile.tracksNA;
            curTracksNA = curTracksNA(idGroupSelectedTogether);
            idGroupSelected = idGroupSelectedReordered;
        catch
            curImportFilePathTracks = fullfile(PathName,'tracksNA.mat');
            curTracksNAfile = load(curImportFilePathTracks,'tracksNA');
            curTracksNA = curTracksNAfile.tracksNA;
            curTracksNA = curTracksNA(idGroupSelectedTogether);
            idGroupSelected = idGroupSelectedReordered;
        end
        mdPath = fileparts(fileparts(fileparts(fileparts(PathName))));
        try
            curMDFile =  load([mdPath filesep 'movieData.mat'],'MD');
        catch
            mdPath = fileparts(fileparts(fileparts(PathName)));
            curMDFile =  load([mdPath filesep 'movieData.mat'],'MD');
        end
        MD = curMDFile.MD;
    end            
end

