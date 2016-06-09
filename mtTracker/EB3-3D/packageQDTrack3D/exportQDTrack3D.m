function exportQDTrack3D(versionPath)

%% Copy code
mkdir([versionPath filesep 'code']);
funList=getFunDependencies('QDTrack3DPackage',{'@Tracks','yalmip','load.m'},true);
for i = 1:numel(funList)
    [p,n,e]=fileparts(funList{i});
    copyfile(funList{i}, [versionPath filesep 'code' filesep n e]);
end

tracksFolder=fileparts(which('Tracks'));
copyfile(tracksFolder,[versionPath filesep 'code' filesep '@Tracks']);
tracksFolder=fileparts(which('kalmanInitLinearMotion'));
copyfile(tracksFolder,[versionPath filesep 'code' filesep 'trackWithGapClosing']);
tracksFolder=fileparts(which('asymDeterm2D3D'));
copyfile(tracksFolder,[versionPath filesep 'code' filesep 'asymDeterm2D3D']);
tracksFolder=fileparts(which('TracksHandle'));
copyfile(tracksFolder,[versionPath filesep 'code' filesep '@TracksHandle']);
tracksFolder=fileparts(which('TracksStruct'));
copyfile(tracksFolder,[versionPath filesep 'code' filesep '@TracksStruct']);
tracksFolder=fileparts(which('lap.m'));
copyfile(tracksFolder,[versionPath filesep 'code' filesep 'linearAssignment']);

copyfile(which('normalizeTracks'),[versionPath filesep 'code']);

%% Copy script and README
userScriptDir=[versionPath filesep 'user-script' filesep];
mkdir(userScriptDir);
movefile([versionPath filesep 'code' filesep 'scriptQuantumDotsTracking.m'],userScriptDir);
movefile([versionPath filesep 'code' filesep 'createQDMovieManagementFile.m'],userScriptDir);
movefile([versionPath filesep 'code' filesep 'simulateAndTrack.m'],userScriptDir);
