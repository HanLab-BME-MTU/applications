function renameSlists
% adds ".mat" to slists

% find path/data
topPath = uigetdir;

% find slists without .mat
slistList = searchFiles('slist_','.mat',topPath,1);

% loop and reassign
for i=1:length(slistList)
    oldFile = [slistList{i,2},filesep,slistList{i,1}];
    movefile(oldFile,[oldFile,'.mat']);
end