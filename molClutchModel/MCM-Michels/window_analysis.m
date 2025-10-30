clc
if usejava('desktop')
    [fileSFolders, pathSFolders] = uigetfile('*.mat','Select MovieData file');
else
    disp({'Enter Path to MovieData (.mat)';
        ['Your current path: ' pwd]});
    rawPath = input(': ','s');
    if isempty(rawPath)
        pathSFolders = 0;
    else
        [pathSFolders, fileSFolders] = fileparts(rawPath);
    end
end

try 
    MD=MovieData.load(append(pathSFolders,fileSFolders));
catch
    disp(['Error :: failed to load file '  fileSFolders])
end

disp(['Loaded MovieData ::' MD.movieDataFileName_])

iFlow = MD.getProcessIndex('FlowAnalysisProcess');
if isempty(iFlow)
    iFlow = MD.getProcessIndex('FlowTrackingProcess');
    if isempty(iFlow)
        error('Flow tracking has to be run.')
    else
        disp('Flow analysis process has not been run, flow tracking process is used instead.')
    end
end
disp(['Analysing...'])

ti=strrep(fileSFolders,'_','\_');
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
%set(t1,'Rotation',90);
%set(t2,'Rotation',90);

hold off


fs=1/ts;
figure()
title(['Power Spectrum of Segment ' num2str(segmentNum)]);
subtitle([ti]);
xlabel("Frequency (Hz)");
ylabel("Power")
hold on
l=length(speed(segmentNum,:));
ff=fft(speed(segmentNum,:));
p=abs(ff).^2;
p1=p(1:floor((n+1)/2));
p1(1)=0;
f=fs/l*(0:floor((l-1)/2));
plot(f,p1);
hold off

runAll=[1:n];
runMaxima=[indexMax,indexMin,indexMaxNext];
selected=[7,9,10,13,15,18,20,22,25];

for i=selected
    s=win.loadChannelOutput(i);
    seg(i)=s(segmentNum);
    seg(i)=seg{i}(1);
    poly{i}=polyshape(cell2mat(seg{i})');
    [cx,cy]=centroid(poly{i});
    figure()
    imshow(cell2mat(images(i)),[125,300])
    xlim([cx-windowSize,cx+windowSize]);
    ylim([cy-windowSize,cy+windowSize]);
    hold on 
    plot(poly{i});
    vectorIn{i}=isinterior(poly{i},flow{i}(:,1),flow{i}(:,2));
    for j=1:length(flow{i})
        for k=1:length(flow{i}(1,:))
        flow{i}(j,k)=flow{i}(j,k)*vectorIn{i}(j);
        end
    end
    quiver(flow{i}(:,1),flow{i}(:,2),vectorScale*flow{i}(:,3),vectorScale*flow{i}(:,4),'AutoScale','off','LineWidth',1,'Color','red')

    sx=(cx-windowSize)+windowSize*0.2;
    sy=(cy+windowSize)-windowSize*0.2;
    quiver(sx,sy,scale10,0,'AutoScale','off','LineWidth',1,'Color','white','MaxHeadSize',5)
    text(sx+scale10*1.2,sy,'500 nm/min','Color','white','FontSize',14)

    fx=(cx-windowSize);
    fy=(cy-windowSize)+windowSize*0.1;
    text(fx+5,fy+5,['Frame ' num2str(i)],'Color','white','FontSize',14)
%     
%     xbar=[sx:sx+scale10];
%     ybar=sy*ones(length(xbar));
%     line(xbar,ybar,'Color','white','LineStyle','-','LineWidth',5)
%     t=text(xbar(end),ybar(1)-1,'100 \mum/min','Color','white');
%     set(t,'Rotation',60);
%     xbar2=[sx:sx+scale10*.75];
%     ybar2=sy*ones(length(xbar2));
%     line(xbar2,ybar2,'Color','blue','LineStyle','-','LineWidth',5)
%     xbar3=[sx:sx+scale10/2];
%     ybar3=sy*ones(length(xbar3));
%     line(xbar3,ybar3,'Color','white','LineStyle','-','LineWidth',5)
%     t2=text(xbar3(end),ybar3(1)-1,'50 \mum/min','Color','white');
%     set(t2,'Rotation',60);
%     xbar4=[sx:sx+scale10/4];
%     ybar4=sy*ones(length(xbar4));
%     line(xbar4,ybar4,'Color','blue','LineStyle','-','LineWidth',5)
%     t3=text(xbar4(1),ybar3(1)-1,'0 \mum/min','Color','white');
%     set(t3,'Rotation',60);


    hold off
end

