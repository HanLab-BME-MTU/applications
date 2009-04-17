function [fullSubIdx]=plusTipPlotRandTraj(projData,type,num,beforeAfterRange,dispType)

% load projData
if nargin<1 || isempty(projData)
    [fileName,pathName]=uigetfile('*.mat','Please select projData from META directory');
    if isequal(fileName,0)
        return
    end
    projData=load([pathName filesep fileName]);
    projData=projData.projData;
end
% set default track type to plot
if nargin<2 || isempty(type)
    type=3;
end
% set default number of tracks to plot
if nargin<3 || isempty(num)
    num=10;
end
% set default range of subtracks before and after the selected subtrack
if nargin<4 || isempty(beforeAfterRange)
    beforeAfterRange=[-inf inf];
end
b=beforeAfterRange(1);
a=beforeAfterRange(2);
% set default for numbered subtracks, legend, or both
if nargin<4 || isempty(dispType)
    dispType=1;
end





% get subtrack coordinates in nSubtracks x nFrames matrices
[xMat,yMat]=plusTipGetSubtrackCoords(projData);

% extract track data
trackData=projData.nTrack_start_end_velMicPerMin_class_lifetime;

% track types: 1 growth, 2 pause/outOfFocus, 3 shrinkage, 4 unclassified
trackType=trackData(:,5);

% displacement per frame (in pixels) matrices
vMat=sqrt(diff(xMat,1,2).^2+diff(yMat,1,2).^2);
vMat(trackType==3,:)=-vMat(trackType==3,:); % make shrink velocity negative

% get indices for num random subtracks of specified type
idx=randsample(find(trackType==type),num);

% full track numbers for all subtracks
fullTrackIdx=trackData(:,1);
% all subtracks in the full track containing idx subtracks
fIdx=arrayfun(@(i) find(fullTrackIdx==fullTrackIdx(i)),idx,'UniformOutput',0);
% first subtrack in full track for each idx subtrack
sIdx=cellfun(@(i) i(1),fIdx);
% last subtrack in full track for each idx subtrack
eIdx=cellfun(@(i) i(end),fIdx);

% modify start/end so they correspond to before/after range
b=max([repmat(b,[num,1]) sIdx-idx],[],2);
a=min([repmat(a,[num,1]) eIdx-idx],[],2);

% combine subtracks into a single row
t=arrayfun(@(i) nansum(vMat(idx(i)+b(i):idx(i)+a(i),:)),[1:num]','UniformOutput',0);
% make nFrames x nTracks matrix
t=cell2mat(t)';
% take cumulative sum to get trajectory
p=cumsum(t);
% before/after track begins, fill with nan
p(t==0)=nan;

% frames where idx tracks begin (first,f) and end (last,l)
f=cell2mat(arrayfun(@(i) find(~isnan(vMat(i,:)),1,'first'),idx,'UniformOutput',0));
l=cell2mat(arrayfun(@(i) find(~isnan(vMat(i,:)),1,'last'),idx,'UniformOutput',0));

% frames where full tracks end
e=cell2mat(arrayfun(@(i) find(~isnan(p(:,i)),1,'last'),[1:num]','UniformOutput',0));

% cumulative displacement value at those frames
vf=p(sub2ind(size(p),f,[1:num]'));
vl=p(sub2ind(size(p),l,[1:num]'));
ve=p(sub2ind(size(p),e,[1:num]'));



colorMap=varycolor(num);
figure
set(gcf,'DefaultAxesColorOrder',colorMap)
h=plot(p);
hold on


% get x,y,color for all points along each track
nFr=size(p,1);
pTemp=p(:);
nanIdx=isnan(pTemp);
allFrm=repmat([1:nFr]',[1 num]);
allFrm=allFrm(:);
allFrm(nanIdx)=[];
pTemp(nanIdx)=[];
col=cell2mat(arrayfun(@(i) repmat(colorMap(i,:),[nFr 1]), [1:num]','UniformOutput',0));
col(nanIdx,:)=[];
% plot all points in p
scatter(allFrm,pTemp,'Marker','.','cData',col)

% plot first and last frames of idx subtracks bigger
scatter(f,vf,'Marker','.','cData',colorMap,'sizedata',round((72/4)^2));
scatter(l,vl,'Marker','.','cData',colorMap,'sizedata',round((72/4)^2));


switch dispType
    case 1
        for i=1:num
            text(e(i),ve(i),['\leftarrow' num2str(fullTrackIdx(idx(i)))],'Color',colorMap(i,:),'FontWeight','bold')
        end
    case 2
        legendEntries=arrayfun(@(i) num2str(fullTrackIdx(idx(i))), [1:num]','UniformOutput',0);
        legend(legendEntries,'location','Best')

    case 3
        for i=1:num
            text(e(i),ve(i),['\leftarrow' num2str(fullTrackIdx(idx(i)))],'Color',colorMap(i,:),'FontWeight','bold')
        end
        legendEntries=arrayfun(@(i) num2str(fullTrackIdx(idx(i))), [1:num]','UniformOutput',0);
        legend(legendEntries,'location','Best')
end

switch type
    case 1
        str='growth';
    case 2
        str='pause';        
    case 3
        str='shrinkage';
    case 4
        str='unclassified';
end
title([num2str(num) ' random ' str ' events, [i-' num2str(abs(beforeAfterRange(1))) ':i+' num2str(abs(beforeAfterRange(2))) ']'])

fullSubIdx=[fullTrackIdx(idx) idx];

