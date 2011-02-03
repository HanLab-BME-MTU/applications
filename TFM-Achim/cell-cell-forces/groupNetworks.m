if isempty(groupedClusters)
    groupedClusters.numClusters= 0;
    groupedClusters.clusterList=[];
end

fileList=searchFiles('trackedNet.mat',[],pwd,1,[],1);
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