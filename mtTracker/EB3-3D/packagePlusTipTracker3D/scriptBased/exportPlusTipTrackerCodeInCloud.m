function exportPlusTipTrackerCodeInCloud(versionNumber)

versionNumberStr=num2str(versionNumber,'%0.2f');
codePath=['/work/gdanuser/proudot/project/EB3-3D-track/packaging/alpha/plusTipTracker3D-alpha-' versionNumberStr];
[p,n,e]=fileparts(codePath);
exportPlusTipTrackerCode(codePath);
system(['cd ' p '; zip -r ' n e '.zip ' n e])
copyfile([codePath '.zip'],[fileparts(codePath) filesep 'linkToCloud' filesep]);