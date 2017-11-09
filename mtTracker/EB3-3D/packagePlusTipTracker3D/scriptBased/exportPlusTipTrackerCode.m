function exportPlusTipTrackerCode(versionPath)

%% Copy code
mkdir([versionPath filesep 'code']);
funList=getFunDependencies('fourScript',{'@Tracks','yalmip','load.m'},true);
for i = 1:numel(funList)
    [p,n,e]=fileparts(funList{i});
    copyfile(funList{i}, [versionPath filesep 'code' filesep n e]);
end

tracksFolder=fileparts(which('Tracks'));
copyfile(tracksFolder,[versionPath filesep 'code' filesep '@Tracks']);
% tracksFolder=fileparts(which('kalmanInitLinearMotion'));
% copyfile(tracksFolder,[versionPath filesep 'code' filesep 'trackWithGapClosing']);
tracksFolder=fileparts(which('TracksHandle'));
copyfile(tracksFolder,[versionPath filesep 'code' filesep '@TracksHandle']);
tracksFolder=fileparts(which('TracksStruct'));
copyfile(tracksFolder,[versionPath filesep 'code' filesep '@TracksStruct']);

copyfile(which('normalizeTracks'),[versionPath filesep 'code']);
copyfile(which('getFeatFromIdx'),[versionPath filesep 'code'])

%% Copy script and README
userScriptDir=[versionPath filesep 'user-script' filesep];
mkdir(userScriptDir);
movefile([versionPath filesep 'code' filesep 'scriptAnalyseResults.m'],userScriptDir);
movefile([versionPath filesep 'code' filesep 'scriptEB3DetectAndTracking.m'],userScriptDir);
movefile([versionPath filesep 'code' filesep 'scriptKinetochoreTracking.m'],userScriptDir);
movefile([versionPath filesep 'code' filesep 'createMovieManagementFile.m'],userScriptDir);
movefile([versionPath filesep 'code' filesep 'loadMoviesManagementFile.m'],userScriptDir);
movefile([versionPath filesep 'code' filesep 'scriptDeskewLatticeData.m'],userScriptDir);
movefile([versionPath filesep 'code' filesep 'exampleCellSelection.m'],userScriptDir);
movefile([versionPath filesep 'code' filesep 'scriptCroppingMovie.m'],userScriptDir);

readmePath=(which('README.txt'));
copyfile(readmePath,versionPath);