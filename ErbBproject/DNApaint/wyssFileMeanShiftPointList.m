function success = wyssFileMeanShiftPointList(PointList,name)
% takes in a directory that contains .mat files with data from the wyss
% converts the data into a PointList structure, drift corrects and aligns
% the data and then runs meanshift clustering on the data.
%
% The .mat files must contain a matrix trace and a matrix drift
%
% Results are then saved in the parent directory of dir as "dir".mat
%

%approx parameters for wyss microscope
Imsize = [256,256];
difLim = 2.1075;

success = 1;

TotalPnts = [];

for i=1:numel(list);
   TotalPnts = vertcat(TotalPnts,[PointList{i}.pnts,i*ones(size(PointList{i}.pnts(:,1)))]); 
end



% since pixelSize of wyss system is ~107 nm this is 32 nm bandwitdh
[clusterInfo,clusterMap]=MeanShiftClustering(TotalPnts(:,1:2),0.3,'kernel','flat','flagDebug',true);

clusterInfo = addToClusterInfo(clusterInfo,TotalPnts,'pixelSize',107);

save([name,'_MeanShiftCluster.mat'],'clusterInfo','clusterMap','TotalPnts');

end

