
[speedCell,protSpeedCell] = quantifyMovieFlowSpeed(MD);
speed=squeeze(speedCell(:,1,:));
[m, n]=size(speed);


flowProcess = MD.getProcess(iFlow);
iChan = find(flowProcess.checkChannelOutput);
flow=flowProcess.loadChannelOutput(iChan,'output','Md');
cf=@(x)[x(:,[2,1]) x(:,[4,3])-x(:,[2,1])];
for i=1:n
    flow{i}=cf(flow{i});
end




iWinPack = MD.getPackageIndex('WindowingPackage');
winPack = MD.getPackage(iWinPack);
win=winPack.getProcess(2);
% use win.loadChannelOutput(t) where t is frame


bf=bfopen(MD.channels_.channelPath_);
images=bf{1}(:);
images=images(1:(length(images)/2));
img=cell(1,length(images));
%imshow(cell2mat(images(1)),[]);

segmentNum=30; %199
windowSize=50;


ts=6;%6 second time step
meta=bf{2};
px=meta.get('spatial-calibration-x'); %um
py=meta.get('spatial-calibration-y');

vectorScale=25;

scale10=vectorScale*str2double(px)*0.5/(ts/60); %500 nm/min

seg=cell(1,n);
poly=cell(1,n);
vectorIn=cell(1,n);

ti=strrep(fileSFolders,'_','\_');
figure()
title({['Speed of Segment ' num2str(segmentNum)]; ti})
% subtitle([ti]);
xlabel('Time (seconds)');
hold on
plot([1:length(speed(segmentNum,:))],speed(segmentNum,:),'LineWidth',2);
s=speed(segmentNum,:);
sTruncMax=s(2:end);
maxs=islocalmax(sTruncMax);
maxsV=sTruncMax(maxs);
indexMax=find(s==maxsV(1));
indexMaxNext=find(s==maxsV(2));
% if rem(length(maxsV),2)==0
%     maxsV=cat(2,maxsV,0);
% end
% medianMax=median(maxsV);
%indexMax=find(s==medianMax);
sTrunc=s(indexMax:end);
mins=islocalmin(sTrunc);
minsV=sTrunc(mins);
indexMin=find(s==minsV(1));
scatter([indexMax,indexMin,indexMaxNext],[s(indexMax),s(indexMin),s(indexMaxNext)],'x','LineWidth',4)
t1=text(indexMax,s(indexMax),['  Frame ' num2str(indexMax)]);
t2=text(indexMin,s(indexMin),['  Frame ' num2str(indexMin)]);
t3=text(indexMaxNext,s(indexMaxNext),['  Frame ' num2str(indexMaxNext)]);