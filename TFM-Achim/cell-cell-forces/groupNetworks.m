function groupedNetworks=groupNetworks(groupedNetworks,filename,doSave,corrTpts)

if nargin<1 || isempty(groupedNetworks)
    groupedNetworks.numClusters= 0;
    groupedNetworks.clusterList=[];
end

if nargin<2 || isempty(filename)
    filename='trackedNet.mat';
end

if nargin<3 || isempty(doSave)
    doSave=1;
end

if nargin<4 || isempty(corrTpts)
    corrTpts=1;
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

rootDir=pwd;

for entryId=1:numel(fileList)
    % load the file
    filestruc=load(fileList{entryId});
    if strcmpi(filename,'trackedNet.mat');
        currNet=filestruc.trackedNet;
    else
        currNet=filestruc.trackedNetCorrected;
    end
    fnameFirstBeadImg=filestruc.fnameFirstBeadImg;
    
    % lines form here...
    if corrTpts
        shrtPathNet=fileList{entryId};
        shrtPathNet=shrtPathNet(2:end);
        fullPathNet=[rootDir,shrtPathNet];
        [pathMechTFM]=getFilenameBody(fullPathNet);
        [projPath]=getFilenameBody(pathMechTFM);
        beadsPath=[projPath,filesep,'data/Beads'];
        beadsFileList=getFileListFromFolder(beadsPath);
        [timePtsRel timePts timeIntervals meanDT stdDT]=getTimeList(beadsFileList,[],1);
        for iframe=1:length(currNet)
            if ~isempty(currNet{iframe})
                currNet{iframe}.par.t      = timePtsRel(iframe);
                currNet{iframe}.par.dt_mean= meanDT;
                currNet{iframe}.par.dt_std = stdDT;
            end
        end        
        if strcmpi(filename,'trackedNet.mat');
            trackedNet=currNet;
            movefile(fullPathNet,[pathMechTFM,filesep,'trackedNetOldDT.mat']);
            save(fullPathNet, 'trackedNet','fnameFirstBeadImg','-v7.3');
        else
            trackedNetCorrected=currNet;
            movefile(fullPathNet,[pathMechTFM,filesep,'trackedNetCorrectedOldDT.mat']);
            save(fullPathNet, 'trackedNetCorrected','fnameFirstBeadImg','-v7.3');
        end
       
    end
    % ... till here could be removed, once time measures in all
    % trackedNet.mat have been corrected.
    
    %save([pwd,filesep,'test.mat'], 'currNet','fnameFirstBeadImg','-v7.3');
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

if doSave && strcmpi(filename,'trackedNet.mat');
    save('groupedNetwork.mat','groupedNetworks');%,'-v7.3')
elseif doSave
    groupedNetworksCorrected=groupedNetworks;
    save('groupedNetworkCorrected.mat','groupedNetworksCorrected');%,'-v7.3')
end