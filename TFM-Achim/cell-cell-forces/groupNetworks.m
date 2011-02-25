function groupedNetworks=groupNetworks(groupedNetworks,filename,doSave)

if nargin<1 || isempty(groupedNetworks)
    groupedNetworks.numClusters= 0;
    groupedNetworks.clusterList=[];
end

if nargin<2 || isempty(filename)
    filename='trackedNet.mat';
end

if isunix
    tic;
    display('Used Unix command find!')
    [~, unixList]=system(['find . -name ''',filename,'''']);
    cellList = textscan(unixList,'%s');
    fileList = cellList{1};
    toc;
else
    % This is another option to search for the files:
    tic;
    display('Couldn''t use unix command find!')
    fileList=searchFiles(filename,[],pwd,1,[],1);
    toc;
end

if nargin<3 || isempty(doSave)
    doSave=1;
end



for entryId=1:numel(fileList)
    % load the file
    filestruc=load(fileList{entryId});
    currNet=filestruc.trackedNet;
    fnameFirstBeadImg=filestruc.fnameFirstBeadImg;
    
    %save([pwd,filesep,'test.mat'], 'currNet','fnameFirstBeadImg','-v7.3');
    
    % sort it in, if it doesn't exist yet
    if groupedNetworks.numClusters>0 && sum(strcmp(fnameFirstBeadImg,groupedNetworks.clusterList))>0
        display('!!! Cluster is already in the group list, nothing to do: ')
        display(fileList{entryId});
        return;
    else
        
        clusterId = groupedNetworks.numClusters+1;
        groupedNetworks.clusterList{clusterId}        = fnameFirstBeadImg;
        groupedNetworks.cluster{clusterId}.trackedNet = currNet;
        groupedNetworks.numClusters                   = groupedNetworks.numClusters+1;
        
        display(['Added: ',fnameFirstBeadImg])
    end
end

if doSave
    save('groupedNetwork.mat','groupedNetworks')
end
