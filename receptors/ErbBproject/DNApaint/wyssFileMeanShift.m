function success = wyssFileMeanShift(dir)
% takes in a directory that contains .mat files with data from the wyss
% converts the data into a PointList structure and then runs meanshift
% clustering on the data.
%
% Results are then saved in the parent directory of dir as "dir".mat
%

success = 1;
cd(dir);
if strcmp(dir(end),filesep)
    dir = dir(1:end-1);
end
output = [dir,'.mat'];

list = what;
list = list.mat;


PointList = [];
TotalPnts = [];

for i = 1:numel(list);

load(list{i});
PointList = vertcat(PointList, {struct('pnts',trace(:,2:3),'name',list{i},'fullData',trace,'drift',[],'dmark',[],'shift',[],'com',[])});
TotalPnts = vertcat(TotalPnts,[trace(:,2:3),i*ones(size(trace(:,1)))]);

end

save([dir,'_PointList.m'],'PointList');

% since pixelSize of wyss system is ~107 nm this is 32 nm bandwitdh
[clusterInfo,clusterMap]=MeanShiftClustering(TotalPnts(:,1:2),0.3,'kernel','flat');

clusterInfo = addToClusterInfo(clusterInfo,TotalPnts,'pixelSize',107);

save([dir,'_MeanShiftCluster.mat'],'clusterInfo','clusterMap','TotalPnts');

end
