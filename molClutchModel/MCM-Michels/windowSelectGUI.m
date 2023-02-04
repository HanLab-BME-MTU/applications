
figure("Name","Window Select","Units","normalized","Position",[.25 .25 .5 .5]);
backPanel=uipanel('BackgroundColor','w');
ws=uipanel('Parent',backPanel,"Title","File Select","FontSize",10,"Units","pixels",'Position',[5,5,150,250]);
ax=axes('Parent',backPanel,'Position',[0.3,0.25,0.75,0.75],"Colormap",gray());
setappdata(backPanel,'img',ax);
set(ax,'XColor','none','YColor','none')
selectedFile=uicontrol(ws,"Style","edit","String","No File Selected","Position",[0 200 150 20],"Enable","inactive","Max",100);
setappdata(backPanel,'file',selectedFile);




movies=uicontrol(ws,"Style","listbox","Position",[0,40,150,120]);
setappdata(backPanel,'selectImg',movies)



fileSelectButton=uicontrol(ws,"String","Select File","Callback",@(btn,event) fileSelectCallback(btn,event),"Position",[0,180,150,20]);

listSelectButton=uicontrol(ws,"String","Load Image","Callback",@(btn,event) loadImage(btn,event),"Position",[0,20,150,20]);


previewGraphs=uipanel('Parent',backPanel,'Position',[0.5,0,0.5,0.25]);
graph1=axes('Parent',previewGraphs,'Position',[0.05,0.1,0.4,0.8]);
graph0=axes('Parent',previewGraphs,'Position',[0.55,0.1,0.4,0.8]);
setappdata(backPanel,'graphs',previewGraphs);


selectedPanel=uipanel('Parent',backPanel,"Title","Selected Windows","FontSize",10,"Units","pixels",'Position',[155,5,150,250]);
windowSelect=uicontrol(selectedPanel,"Style","listbox","Position",[0,40,150,120],"Enable","inactive");
setappdata(backPanel,'windows',windowSelect);
saveSelectedButton=uicontrol(selectedPanel,"String","Save Selected","Callback",@(btn,event) saveSelectedCallback(btn,event),"Position",[0,20,150,20]);
loadSelectButton=uicontrol(selectedPanel,"String","Select Data File","Callback",@(btn,event) dataSelectCallback(btn,event),"Position",[0,180,150,20]);
selectedData=uicontrol(selectedPanel,"Style","edit","String","No Data Selected","Position",[0 200 150 20],"Enable","inactive","Max",100);
setappdata(backPanel,'selectedData',selectedData);

sliderPanel=uipanel('Parent',backPanel,'Position',[0.3,0.25,0.2,0.05]);
maxSlider=uicontrol('Parent',sliderPanel,'Style','slider','Units','normalized','Position',[0,.5,1,.5],'Value',.5,"Callback",@(src,event) maxImageCallback(src,event));
minSlider=uicontrol('Parent',sliderPanel,'Style','slider','Units','normalized','Position',[0,0,1,.5],'Value',.1,"Callback",@(src,event) minImageCallback(src,event));

data=struct;
setappdata(backPanel,'data',data);






function fileSelectCallback(btn,ax)
    sf=getappdata(btn.Parent.Parent,'file');
    fig=getappdata(btn.Parent.Parent,"img");
    list=getappdata(btn.Parent.Parent,"selectImg");
    data=getappdata(btn.Parent.Parent,"data");
    [fileSFolders, pathSFolders] = uigetfile('*.mat','Select MovieList file');

    sf.String=fileSFolders;
    setappdata(sf,'path',append(pathSFolders,fileSFolders))

    try 
        ML=load(append(pathSFolders,fileSFolders));
    catch
        disp(['Error :: failed to load file '  fileSFolders])
    end
    ML=ML.ML;
    disp(['Loaded MovieList ::' ML.movieListFileName_])
    list.String=ML.movieDataFile_;
    for i = [1:length(ML.movieDataFile_)]
        data(i).file=ML.movieDataFile_{i};
        data(i).selected={};
    end
    setappdata(btn.Parent.Parent,'data',data);
    setappdata(btn.Parent.Parent,'movieList',ML.movieDataFile_);
end
function loadImage(btn,ax)
    ml=getappdata(btn.Parent.Parent,'movieList');
    fig=getappdata(btn.Parent.Parent,"img");
    list=getappdata(btn.Parent.Parent,"selectImg");
    windowList=getappdata(btn.Parent.Parent,"windows");
    data=getappdata(btn.Parent.Parent,"data");

       

    l=cellfun(@num2str,data(list.Value).selected,'UniformOutput',false);
    windowList.String=l;
    try 
        MD=MovieData.load(ml{list.Value});
    catch
        disp(['Error :: failed to load file '  ml{list.Value}])
    end
    
    disp(['Loaded MovieData ::' MD.movieDataFileName_])
    disp(['Analysing...'])
    
    [speedCell,protSpeedCell] = quantifyMovieFlowSpeed(MD);
    speed=squeeze(speedCell(:,1,:));
    [m, n]=size(speed);

    iWinPack = MD.getPackageIndex('WindowingPackage');
    winPack = MD.getPackage(iWinPack);
    win=winPack.getProcess(2);
    bf=bfopen(MD.channels_.channelPath_);
    images=bf{1}(:);
    images=images(1:(length(images)/2));

    frame=1;
    image(cell2mat(images(frame)),'Parent',fig,'CDataMapping','scaled');
    fig.Colormap=gray();
    fig.CLim=[100, 500];
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
            
                sel=logical(sum(cell2mat(cellfun(@(x)isequal(x,[j l]),data(list.Value).selected,'UniformOutput',false))));
                if sel
                    col='g';
                else
                    col='r';
                end
                p=plot(fig,poly{frame},"FaceColor",col,'EdgeColor',col,'ButtonDownFcn',@(src,event)getSelected(src,event));
                ind=[j,l];
                pSpeed=squeeze(speedCell(j,l,:));
                setappdata(p,'index',ind)
                setappdata(p,'speed',pSpeed)
                setappdata(p,'sel',sel)
            end
            %drawnow;
        end

end


function getSelected(src,event)

    fig=getappdata(src.Parent.Parent,"graphs");
    g0=fig.Children(2);
    g1=fig.Children(1);
    ind=getappdata(src,'index');
    speed=getappdata(src,'speed');
    sel=getappdata(src,'sel');
    listFig=getappdata(src.Parent.Parent,"selectImg");
    list=listFig.Value;
    data=getappdata(src.Parent.Parent,"data");
    windowList=getappdata(src.Parent.Parent,"windows");
    %figure("Name","Speed and Spectrum","Units","normalized","Position",[0.25,0.25,0.25,0.5]);
    l=length(speed);
    fs=(1/6);
    ff=fft(speed-mean(speed));
    p=abs(ff).^2;
    p1=p(1:floor((l+1)/2));
    f=fs/l*(0:floor((l-1)/2));
    
    plot(g0,[1:length(speed)]*6,speed);
    
    g0.Title.String=(['Speed (Window : [',num2str(ind(1)),',',num2str(ind(2)),'])']);
    g0.XLabel.String=('Time (seconds)');
    g0.YLabel.String=('Speed');

    plot(g1,f,p1);
    g1.Title.String=(['Power Spectrum of Speed (Window : [',num2str(ind(1)),',',num2str(ind(2)),'])']);
    g1.XLabel.String=('Power');
    g1.YLabel.String=('Frequency');
    if ~sel
        src.FaceColor=[0 1 0];
        src.EdgeColor=[0 1 0];
        setappdata(src,'sel',true);
        data(list).selected{end+1}=[ind(1) ind(2)];
        
    else
        src.FaceColor=[1 0 0];
        src.EdgeColor=[1 0 0];
        setappdata(src,'sel',false)
        data(list).selected(cell2mat(cellfun(@(x)isequal(x,[ind(1) ind(2)]),data(list).selected,'UniformOutput',false)))=[];
    end
    l=cellfun(@num2str,data(list).selected,'UniformOutput',false);
    windowList.String=l;
    setappdata(src.Parent.Parent,'data',data);
end

function saveSelectedCallback(btn,event)
    sf=getappdata(btn.Parent.Parent,'file');
    fp=getappdata(sf,'path');
    data=getappdata(btn.Parent.Parent,"data");
    [f,p]=uiputfile({'*.mat'},'Save as','selected.mat');
    save(append(p,f),'data',"fp")

end



function dataSelectCallback(btn,ax)
    sf=getappdata(btn.Parent.Parent,'selectedData');
    sf1=getappdata(btn.Parent.Parent,'file');
    fig=getappdata(btn.Parent.Parent,"img");
    list=getappdata(btn.Parent.Parent,"selectImg");
    data=getappdata(btn.Parent.Parent,"data");
    [fileSFolders, pathSFolders] = uigetfile('*.mat','Select data file');

    sf.String=fileSFolders;
    
    


    try 
        dataFile=load(append(pathSFolders,fileSFolders));
    catch
        disp(['Error :: failed to load file '  fileSFolders])
    end
    data=dataFile.data;

    try 
        ML=load(dataFile.fp);
    catch
        disp(['Error :: failed to load file '  fileSFolders])
    end
    ML=ML.ML;
    disp(['Loaded MovieList ::' ML.movieListFileName_])
    list.String=ML.movieDataFile_;
    setappdata(btn.Parent.Parent,'movieList',ML.movieDataFile_);
    


    setappdata(sf,'path',dataFile.fp)
    setappdata(btn.Parent.Parent,'data',data);


end



function maxImageCallback(src,event)
    img=getappdata(src.Parent.Parent,'img');
    percent=get(src,'Value');
    img.CLim=[img.CLim(1),percent*1000];
    setappdata(src.Parent.Parent,'img',img)

    
end


function minImageCallback(src,event)
    img=getappdata(src.Parent.Parent,'img');
    percent=get(src,'Value');
    img.CLim=[percent*1000,img.CLim(2)];
    setappdata(src.Parent.Parent,'img',img)
    
end