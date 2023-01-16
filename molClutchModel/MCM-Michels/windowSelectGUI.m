
figure("Name","Window Select");
backPanel=uipanel('BackgroundColor','w');
ws=uipanel('Parent',backPanel,"Title","File Select","FontSize",10,'Position',[0,0.7,0.3,0.2]);
ax=axes('Parent',backPanel,'Position',[0.3,0.25,0.75,0.75]);
setappdata(backPanel,'img',ax);
set(ax,'XColor','none','YColor','none')
selectedFile=uicontrol(ws,"Style","edit","String","No File Selected","Position",[0 40 150 20],"Enable","inactive","Max",100);
setappdata(backPanel,'file',selectedFile);
fileSelectButton=uicontrol(ws,"String","Select File","Callback",@(btn,event) fileSelectCallback(btn,event),"Position",[0,20,150,20]);








function fileSelectCallback(btn,ax)
    sf=getappdata(btn.Parent.Parent,'file');
    fig=getappdata(btn.Parent.Parent,"img");
    [fileSFolders, pathSFolders] = uigetfile('*.mat','Select MovieData file');

    sf.String=fileSFolders;

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

    frame=1;
    image(cell2mat(images(frame)),'Parent',fig);
    set(fig,'XColor','none','YColor','none')
    hold(fig,"on")
    s=win.loadChannelOutput(frame);



    [m, n]=size(speed);
        for j=1:m
            
            seg(j,frame)=s(j);
            for l=1:length(seg{j,frame})
                if l>5
                    continue
                end
                seg(j,frame)=s(j);
                seg(j,frame)=seg{j,frame}(l);
                poly{frame}=polyshape(cell2mat(seg{j,frame})');
                
                
                %xlim([cx-windowSize,cx+windowSize]);
                %ylim([cy-windowSize,cy+windowSize]);
            
                p=plot(fig,poly{frame},"FaceColor",'r','EdgeColor','r','ButtonDownFcn',@(src,event)getSelected(src,event));
                ind=[j,l];
                pSpeed=squeeze(speedCell(j,l,:));
                setappdata(p,'index',ind)
                setappdata(p,'speed',pSpeed)
            end
            drawnow;
        end

end


function getSelected(src,event)
    ind=getappdata(src,'index');
    speed=getappdata(src,'speed');
    disp(ind)
    disp(speed)
    figure("Name","Speed and Spectrum","Units","normalized","Position",[0.25,0.25,0.25,0.5]);
    l=length(speed);
    fs=(1/6);
    ff=fft(speed-mean(speed));
    p=abs(ff).^2;
    p1=p(1:floor((l+1)/2));
    f=fs/l*(0:floor((l-1)/2));


    subplot(1,2,1);
    plot([1:length(speed)]*6,speed);
    title(['Speed (Window : [',num2str(ind(1)),',',num2str(ind(2)),'])'])
    xlabel('Time (seconds)')
    ylabel('Speed')
    subplot(1,2,2);
    plot(f,p1);
    title(['Power Spectrum of Speed (Window : [',num2str(ind(1)),',',num2str(ind(2)),'])'])
    xlabel('Power')
    ylabel('Frequency')
    


end

