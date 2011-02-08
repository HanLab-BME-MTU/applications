function groupedClusters=groupNetworks(groupedClusters,filename)

if nargin<1 || isempty(groupedClusters)
    groupedClusters.numClusters= 0;
    groupedClusters.clusterList=[];
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



for entryId=1:numel(fileList)
    % load the file
    filestruc=load(fileList{entryId});
    currNet=filestruc.trackedNet;
    fnameFirstBeadImg=filestruc.fnameFirstBeadImg;
    
    %save([pwd,filesep,'test.mat'], 'currNet','fnameFirstBeadImg','-v7.3');
    
    % sort it in, if it doesn't exist yet
    if groupedClusters.numClusters>0 && sum(strcmp(fnameFirstBeadImg,groupedClusters.clusterList))>0
        display('Cluster is already in the group list, nothing to do: ')
        display(fileList{entryId});
        return;
    else
        
        clusterId = groupedClusters.numClusters+1;
        groupedClusters.clusterList{clusterId}        = fnameFirstBeadImg;
        groupedClusters.cluster{clusterId}.trackedNet = currNet;
        groupedClusters.numClusters                   = groupedClusters.numClusters+1;
    end
end