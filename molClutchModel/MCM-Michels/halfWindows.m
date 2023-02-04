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
c=jet(5);
layers=5;
i=1;
    figure()
    imshow(cell2mat(images(i)),[125,300])
    hold on 
        s=win.loadChannelOutput(i);
        speed=squeeze(speedCell(:,1,:));
        [m, n]=size(speed);
        for j=1:floor(m/2)
            
            seg(j,i)=s(j);
            for l=1:length(seg{j,i})
                if l>5
                   continue
                end
                seg(j,i)=s(j);
                seg(j,i)=seg{j,i}(l);
                poly{i}=polyshape(cell2mat(seg{j,i})');
                
                
                %xlim([cx-windowSize,cx+windowSize]);
                %ylim([cy-windowSize,cy+windowSize]);
            
                plot(poly{i},"FaceColor",c(6-l,:),'EdgeColor',c(6-l,:));
            end
            drawnow;
        end

%     vectorIn{i}=isinterior(poly{i},flow{i}(:,1),flow{i}(:,2));
%     for j=1:length(flow{i})
%         for k=1:length(flow{i}(1,:))
%         flow{i}(j,k)=flow{i}(j,k)*vectorIn{i}(j);
%         end
%     end
%     quiver(flow{i}(:,1),flow{i}(:,2),vectorScale*flow{i}(:,3),vectorScale*flow{i}(:,4),'AutoScale','off','LineWidth',1,'Color','red')
% 
%     sx=(cx-windowSize)+windowSize*0.2;
%     sy=(cy+windowSize)-windowSize*0.2;
%     quiver(sx,sy,scale10,0,'AutoScale','off','LineWidth',1,'Color','white','MaxHeadSize',5)
%     text(sx+scale10*1.2,sy,'500 nm/min','Color','white','FontSize',14)
% 
%     fx=(cx-windowSize);
%     fy=(cy-windowSize)+windowSize*0.1;
%     text(fx+5,fy+5,['Frame ' num2str(i)],'Color','white','FontSize',14)
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


