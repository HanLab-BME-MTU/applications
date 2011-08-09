function [path]=findPath2FileName(filename)

if isunix
    tic;
    display('Used Unix command find!')
    [~, unixList]=system(['find . -name ''',filename,'''']);
    cellList = textscan(unixList,'%s');
    path = cellList{1};
    toc;
else
    % This is another option to search for the files:
    tic;
    display('Couldn''t use unix command find!')
    path = searchFiles(filename,[],pwd,1,[],1);
    toc;
end

end

% This is to find all used clusters:
% 
% for clusterID=1:length(groupedNetworks.cluster)
%     filename  =groupedNetworks.clusterList{clusterID}
%     [currPath]=findPath2FileName(filename)
%     if length(currPath)>1
%         display(['The path for cluster : ',num2str(clusterID),' is not unique'])
%     end
%     pathList{clusterID}=currPath;
% end