function success = wyssFileMeanShiftPlusDrift(dir)
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
pnts = [pnts,trace(:,10)];
%dmark = findDriftM(pnts(trace(:,1)==1,:),pnts(trace(:,1)==max(trace(:,1)),:),1000);
%dmark = findDriftM2(pnts(trace(:,1)==1,:),pnts,max(trace(:,1)));
dmark = findDriftFit(pnts);
if isempty(dmark)
    error(['missing drift marker ',list{i}]);
end
PointList = vertcat(PointList, {struct('pnts',pnts(:,1:2),'name',list{i},'fullData',trace,'drift',drift,'dmark',dmark,'shift',[],'com',[])});

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
     [shift,transform,A,B] = driftMarkerRegistration(test,ref,Imsize,difLim);
     TotalShift = shift - transform.trans;
     C = test - repmat(TotalShift,[numel(test(:,1)),1]);
     PointList{j}.shift = struct('transform',transform,'preshift',shift,'A',A,'B',B,'Postshift',C,'TotalShift',TotalShift);
     
     %shifts points to align
     tmp = PointList{j}.pnts(:,1:2)-repmat(TotalShift,[numel(PointList{j}.pnts(:,1)),1]);
     PointList{j}.pnts=tmp(~isnan(tmp(:,1)),:);
     
end

 %Since all other Images are shifted to match this one, here only NaNs are
 %removed
 tmp = PointList{j_ref}.pnts(:,1:2);
 PointList{j_ref}.pnts=tmp(~isnan(tmp(:,1)),:);
 
 
save([dir,'_DriftPlus_PointList.mat'],'PointList');

for i=1:numel(list);
   TotalPnts = vertcat(TotalPnts,[PointList{i}.pnts,i*ones(size(PointList{i}.pnts(:,1)))]); 
end



% since pixelSize of wyss system is ~107 nm this is 32 nm bandwitdh
[clusterInfo,clusterMap]=MeanShiftClustering(TotalPnts(:,1:2),0.3,'kernel','flat','flagDebug',true);

clusterInfo = addToClusterInfo(clusterInfo,TotalPnts,'pixelSize',107);

save([dir,'_DriftPlus_MeanShiftCluster.mat'],'clusterInfo','clusterMap','TotalPnts');

end

function dmark=findDriftM(firstF,lastF,minInten)
%identifies drift markers as any points that appear in both the first and
%last frame within two pixels. returns the position in the first frame

firstF(firstF(:,3)<minInten,:)=NaN;
lastF(lastF(:,3)<minInten,:)=NaN;


dm = distMat2(lastF(:,1:2),firstF(:,1:2));
[x,y] = ind2sub(size(dm),find(dm<0.5));
dmark = firstF(y,1:2);

end


function dmark=findDriftM2(firstF,pnts,nFrames)


  img = hist3(pnts(:,1:2),'Edges',{[0:256],[0:256]});

% finds pixels with large localizations nearly equal to the lenght of the movie
% these should be the drift markers

%sets the requirement to be a drift marker to be 200 less than the max value
% if max is less than nFrames-500

maxN = max(img(:));
if maxN < nFrames-500
nFrames = maxN + 300;
end

[x,y] = ind2sub(size(img),find(img>=(nFrames-500)));  

 n= numel(x);
 dmark = zeros(n,2);

% this loop finds the points in the first frame near the identified drift markers
 for i=1:n
		 
		 ind = find(firstF(:,1)>=x-1 & firstF(:,1)<=x+1 & firstF(:,2)>=y-1 & firstF(:,2)<= y+1);
		 dmark(i,:)=mean(firstF(ind,:),1);
 end

end

%yet another function for finding drift markers. This one creates an image via histograming and finds the dmarks by fitting to a 2D gaussian

function [dmark] = findDriftFit(pnts)

   img = hist3(pnts(:,1:2),'Edges',{[0:0.1:256],[0:0.1:256]});

%finds maxium and sets cutoff to 0.95*max
   m =0.9*max(img(:));
   
[x,y]=ind2sub(size(img),find(img >= m));

n = numel(x);

if n >20
  error;
end

  %There is always at least one
dmark = zeros([n,2]);
for i=1:n
	temp = fitGaussians2D(img(x(i)-6:x(i)+6,y(i)-6:y(i)+6),7,7,m,1,0,'xyAcs');        
dmark(i,:)=[(x(i)+(temp.x-7))/10,(y(i)+(temp.y-7))/10];
end

end
