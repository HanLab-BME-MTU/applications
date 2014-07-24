function success = wyssFileMeanShiftManualAligned(PointList, name,varargin)
% takes in a directory that contains .mat files with data from the wyss
% converts the data into a PointList structure, drift corrects and aligns
% the data and then runs meanshift clustering on the data.
%
% PointList is the algined point list of each channel
% name is the name for saving the results
%
% Results are then saved in the parent directory of dir as "dir".mat
%

ip=inputParser;
ip.CaseSensitive=false;
ip.StructExpand=true;

ip.addRequired('PointList',@iscell);
ip.addRequired('name',@ischar);

ip.addOptional('bandW',0.3,@isscalar);

ip.parse(PointList,name,varargin{:});

bandW = ip.Results.bandW;


%approx parameters for wyss microscope
Imsize = [256,256];
difLim = 2.1075;

success = 1;

TotalPnts =[];

for i=1:numel(PointList);
   TotalPnts = vertcat(TotalPnts,[PointList{i}.pnts,i*ones(size(PointList{i}.pnts(:,1)))]); 
end



% since pixelSize of wyss system is ~107 nm this is 32 nm bandwitdh
[clusterInfo,clusterMap]=MeanShiftClustering(TotalPnts(:,1:2),bandW,'kernel','flat','flagDebug',true);

clusterInfo = addToClusterInfo(clusterInfo,TotalPnts,'pixelSize',160.5);

save([name,'_MeanShiftCluster.mat'],'clusterInfo','clusterMap','TotalPnts');

end
