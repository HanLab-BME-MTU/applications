function success = wyssFileMeanShiftTracking(dir, UseMask)
% takes in a directory that contains .mat files with data from the wyss
  % converts the data into a PointList structure, tracks, aligns all channels  and then runs meanshift
% clustering on the data.
% 
%If use Mask is set to 1 then the .mat files must contain a file mask as well as trace
%This mask will identify a region on the movie that contains a drift marker. If set to 1 will use an amplitude threshold to reduce the 
% data size before tracking and then identify long tracks as driftmarkers.
% 
% creates a results directory so that these files can be processed by exchangePaintAlignment.
% 
% Results are then saved in the parent directory of dir as "dir".mat
%

success = 1;
cd(dir);
if strcmp(dir(end),filesep)
    dir = dir(1:end-1);
end
output = [dir,'.mat'];
ind = findstr(dir,filesep)
  name = [dir(ind(end-1):ind(end)-1),'_',dir(ind(end)+1:end)];
list = what;
list = list.mat;


PointList = [];
TotalPnts = [];

for i = 1:numel(list);

load(list{i});
%PointList = vertcat(PointList, {struct('pnts',trace(:,2:3),'name',list{i},'fullData',trace,'drift',[],'dmark',[],'shift',[],'com',[])});

% creates an abbreviated pointList to use by wyssFileMeanShiftTracking
ind = trace(:,10)>5000;
FN = wyssFileTrackingSingle(trace(ind,:),list{i});
load(FN);
%adjust tracksFinal to include each localization as an stand alone "track"

tracksFinal = appendSingleTracks(tracksFinal,trace);
save(FN,'tracksFinal','features');

end

trackingList = findFilesInSubDirs(cd(),'tracking');
fullpath = exchangePaintAlignment(trackingList,[name,'_PointListTracking'],'ImSize',[256,256],'MinTrackLen',1,'MinDriftLen',10);

load(fullpath);

%TotalPnts = vertcat(TotalPnts,[trace(:,2:3),i*ones(size(trace(:,1)))]);

n = numel(PointList);
for i = 1:n
    TotalPnts = vertcat(TotalPnts,[PointList{i}.com,i*ones(size(PointList{i}.com))]);
end

TotalPnts=TotalPnts(~isnan(TotalPnts(:,1)),:);


%save([dir,'_PointList.m'],'PointList');

% since pixelSize of wyss system is ~107 nm this is 32 nm bandwitdh
[clusterInfo,clusterMap]=MeanShiftClustering(TotalPnts(:,1:2),0.3,'kernel','flat');

clusterInfo = addToClusterInfo(clusterInfo,TotalPnts,'pixelSize',107);

save([dir,'_MeanShiftClusterTracking.mat'],'clusterInfo','clusterMap','TotalPnts');

end


function tracksFinal = appendSingleTracks(tracksFinal,trace)
%% Note as is drift markers will still be included in the image.
  numpnts = cellfun(@numel,{tracksFinal.tracksFeatIndxCG});
  ind = numpnts > 1000;
tracksFinal = tracksFinal(ind);

parfor i= 1:numel(trace(:,1));
tracksFinal = vertcat(tracksFinal,struct('tracksFeatIndxCG',i,'tracksCoordAmpCG',[trace(i,2),trace(i,3),0,trace(i,10),trace(i,5),trace(i,6),0,trace(i,9)],'seqOfEvents',[trace(i,1),1,1,NaN;trace(i,1),2,1,NaN]));
end
end
