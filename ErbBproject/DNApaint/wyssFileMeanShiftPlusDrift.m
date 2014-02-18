function success = wyssFileMeanShiftPlusDrift(dir)
% takes in a directory that contains .mat files with data from the wyss
% converts the data into a PointList structure, drift corrects and aligns
% the data and then runs meanshift clustering on the data.
%
% The .mat files must contain a matrix trace and a matrix drift
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
%drift correct using trace
pnts=trace(:,2:3)-drift(trace(:,1),:);
dmark = findDriftM(pnts(trace(:,1)==1,:),pnts(trace(:,1)==max(trace(:,1)),:));
if isempty(dmark)
    error(['missing drift marker ',list{i}]);
end
PointList = vertcat(PointList, {struct('pnts',pnts,'name',list{i},'fullData',trace,'drift',drift,'dmark',dmark,'shift',[],'com',[])});

end

%aligns images a la exchangePaintAlignment (copied from). Only
%translational alignment

k = 1:numel(PointList);
j_ref = 1;
for j=k
    if numel(PointList{j}.dmark(:,1)) > numel(PointList{j_ref}.dmark(:,1))
        j_ref = j;
    end
end

k(j_ref)=[];

ref = PointList{j_ref}.dmark;
s_ref = size(ref);
PointList{j_ref}.shift = struct('transform',[],'preshift',[],'A',[],'B',[],'Postshift',ref);

for j = k
     test = PointList{j}.dmark;
     [shift,transform,A,B] = driftMarkerRegistration(test,ref,ImSize,difLim);
     TotalShift = shift - transform.trans;
     C = test - repmat(TotalShift,[numel(test(:,1)),1]);
     
     %shifts points to align
     tmp = PointList{j}.pnts(:,1:2)-repmat(TotalShift,[numel(PointList{j}.pnts(:,1)),1]);
     PointList{j}.pnts=tmp(~isnan(tmp(:,1)),:);
     
end

 %Since all other Images are shifted to match this one, here only NaNs are
 %removed
 tmp = PointList{j_ref}.pnts(:,1:2);
 PointList{j_ref}.pnts=tmp(~isnan(tmp(:,1)),:);
 
 
save([dir,'_DriftPlus_PointList.m'],'PointList');

for i=1:numel(list);
   TotalPnts = vertcat(TotalPnts,[PointList{i}.pnts,i*ones(size(PointList{i}.pnts(:,1)))]); 
end



% since pixelSize of wyss system is ~107 nm this is 32 nm bandwitdh
[clusterInfo,clusterMap]=MeanShiftClustering(TotalPnts(:,1:2),0.3,'kernel','flat');

clusterInfo = addToClusterInfo(clusterInfo,TotalPnts,'pixelSize',107);

save([dir,'_DriftPlus_MeanShiftCluster.mat'],'clusterInfo','clusterMap','TotalPnts');

end

function dmark=findDriftM(firstF,lastF)
%identifies drift markers as any points that appear in both the first and
%last frame within two pixels. returns the position in the first frame

dm = distMat2(lastF,firstF);
[x,y] = ind2sub(size(dm),find(dm<0.5));
dmark = firstF(y,:);

end

