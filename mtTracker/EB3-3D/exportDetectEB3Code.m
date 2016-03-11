function exportDetectEB3Code(versionPath)

%% Copy code
mkdir([versionPath filesep 'code']);
funList=getFunDependencies('detectPackage',{'@Tracks','yalmip','load.m'},true);
for i = 1:numel(funList)
    [p,n,e]=fileparts(funList{i});
    copyfile(funList{i}, [versionPath filesep 'code' filesep n e]);
end


%% Copy script and README
userScriptDir=[versionPath filesep 'user-script' filesep];
mkdir(userScriptDir);
movefile([versionPath filesep 'code' filesep 'scriptEB3Detect.m'],userScriptDir);
movefile([versionPath filesep 'code' filesep 'createMovieManagementFile.m'],userScriptDir);
movefile([versionPath filesep 'code' filesep 'loadMoviesManagementFile.m'],userScriptDir);

readmePath=(which('README.txt'));
copyfile(readmePath,versionPath);