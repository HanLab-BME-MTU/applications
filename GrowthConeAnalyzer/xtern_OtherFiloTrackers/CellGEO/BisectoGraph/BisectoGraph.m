function varargout = BisectoGraph(varargin)
% BISECTOGRAPH MATLAB code for BisectoGraph.fig
%      BISECTOGRAPH, by itself, creates a new BISECTOGRAPH or raises the existing
%      singleton*.
%
%      H = BISECTOGRAPH returns the handle to a new BISECTOGRAPH or the handle to
%      the existing singleton*.
%
%      BISECTOGRAPH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BISECTOGRAPH.M with the given input arguments.
%
%      BISECTOGRAPH('Property','Value',...) creates a new BISECTOGRAPH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BisectoGraph_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BisectoGraph_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BisectoGraph

% Last Modified by GUIDE v2.5 01-Mar-2013 00:54:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BisectoGraph_OpeningFcn, ...
                   'gui_OutputFcn',  @BisectoGraph_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

function BisectoGraph_OpeningFcn(hObject, eventdata, handles, varargin)

    handles.output = hObject;    
    def=importdata('Defaults.dat');
    addpath([cd '/bg']);

    handles.NUM_ELEMENTS = 0;
    handles.IMAGE_NUMBER = 1;
    handles.X = 0; 
    handles.Y = 0; 
    
    handles.cLINE = 0;
    
    handles.Lag = 1;
    handles.NeedSpeed = 0;
    handles.DefCrR = def(1);
    handles.DefGN = def(2);
    
    handles.CutOff = 0;

    % Listeners
    handles.LisMov = addlistener(handles.movieSlider, 'Value', 'PostSet', @(src, event)switchFrame(hObject, src, event));
    handles.LisLag = addlistener(handles.lagSlider, 'Value', 'PostSet', @(src, event)switchLag(hObject, src, event));
    handles.LisFil = addlistener(handles.CrR_Slider, 'Value', 'PostSet', @(src, event)switchFilt(hObject, src, event));

    set(handles.frame_text,'ForegroundColor',[.5 .5 .5]);
    set(handles.Tip_text,'String','  ');
    set(handles.CrR_text,'String','  ');
    set(handles.results_panel,'ForegroundColor',[.5 .5 .5]);
    set(handles.smoo_text,'ForegroundColor',[.5 .5 .5]);
    set(handles.label_text,'ForegroundColor',[.5 .5 .5]);

    % Update handles structure
    set(handles.boundary_axis,'Box','on','XTick',[],'YTick',[]);
    set(handles.length_width_axis,'Box','on','XTick',[],'YTick',[]);
    set(handles.profile_axis,'Box','on','XTick',[],'YTick',[]);

    set(handles.movieSlider, 'Enable', 'off');
    set(handles.lagSlider, 'Enable', 'off');
    set(handles.CrR_Slider, 'Enable', 'off');
    
    set(handles.graph_frame_button, 'Enable', 'off');
    set(handles.graph_all_button, 'Enable', 'off');
    set(handles.profile_button, 'Enable', 'off');
    set(handles.batch_button, 'Enable', 'off');
    set(handles.list_button, 'Enable', 'off');

    set(handles.checkbox1, 'Enable', 'off');
    
    set(handles.save_res, 'Enable', 'off');
    
    set(handles.uitable1,'Data',cell(1,2),'Enable','off');

    guidata(hObject, handles);

function varargout = BisectoGraph_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function main_window_CloseRequestFcn(hObject, eventdata, handles)

    try 
        def = [handles.CrR, handles.GN];
        id = fopen([cd '/Defaults.dat'],'w');        
        fprintf(id,'%d \n',def);
        fclose(id);
    end
    delete(hObject);    
%--------------------------------------------------------------------------



% --- Function for Listeners ---------------------------------------------- 
function switchFrame(hObject, src, event)

    handles = guidata(hObject);    
    sliderNumber = round(event.NewValue);    
    handles.IMAGE_NUMBER = sliderNumber;
    
    %axes(handles.length_width_axis);
    cla(handles.length_width_axis); 
    %axis fill;    
    set(handles.length_width_axis,'Box','on','XTick',[],'YTick',[],'Visible','on');
        
    %axes(handles.profile_axis);
    cla(handles.profile_axis);
    %axis fill;    
    set(handles.profile_axis,'Box','on','XTick',[],'YTick',[],'Visible','on');
    
    axes(handles.boundary_axis);
    xlim=get(handles.boundary_axis,'XLim');
    ylim=get(handles.boundary_axis,'YLim');
    cla;
    hold on;
    bound=handles.BOUNDARIES{sliderNumber};
    by=bound(:,1); bx=bound(:,2); 
    fill(bx,by,[.9 .9 .9],'EdgeColor','k');
    axis image;
    axis ij;
    set(handles.boundary_axis,'XLim',xlim);
    set(handles.boundary_axis,'YLim',ylim);
    axis off;
    
    
    
    set(handles.frame_text,'String',['Frame # ' num2str(sliderNumber)]);
    
    guidata(hObject, handles);

function switchLag(hObject, src, event) 

    handles = guidata(hObject);
    lag = round(event.NewValue);
    
    old = handles.OLD; 
    if lag~=old   
    
        fr = round(get(handles.movieSlider,'Value'));          
        set(handles.Tip_text, 'String', ['Tip # ' num2str(lag)]);

        if handles.DonePR(fr) 

            Tip=handles.Tip{fr};
            Pr=handles.Pr{fr};
            Py=handles.Py{fr};
            rpath=handles.RedPaths{fr};
            lwcrv=handles.LWcurves{fr}; 
            
            [xbase,ybase,lbase,ibase]=FindBase(rpath{lag},lwcrv{lag},handles.CrR,handles.CutOff);
            
            axes(handles.boundary_axis);
            xlim=get(handles.boundary_axis,'XLim');
            ylim=get(handles.boundary_axis,'YLim');
            cla;
            hold on;
            bound=handles.BOUNDARIES{fr};
            by=bound(:,1); bx=bound(:,2); 
            fill(bx,by,[.9 .9 .9],'EdgeColor','k');
            plot(rpath{lag}(:,1),rpath{lag}(:,2),'r','LineWidth',2);
            plot(bx(Tip),by(Tip),'bo','MarkerFaceColor','b','MarkerSize',2);
            axis image;
            axis ij;
            set(handles.boundary_axis,'XLim',xlim);
            set(handles.boundary_axis,'YLim',ylim);
            axis off;
                        
            axes(handles.profile_axis);
            xlim=get(handles.profile_axis,'XLim');
            ylim=get(handles.profile_axis,'YLim');            
            cla; hold on;            
            plot(Pr,'c');
            plot(Py,'k','LineWidth',1);
            plot(Tip,Py(Tip),'bo','MarkerFaceColor','b','MarkerSize',4);
            set(handles.profile_axis,'Box','on','XLim',xlim,'YLim',ylim); 
            plot([Tip(lag) Tip(lag)],ylim,'r','LineWidth',2);      
            
            axes(handles.length_width_axis);
            ylim=get(handles.length_width_axis,'YLim');
            cla reset; hold on;
            plot(lwcrv{lag}(:,1),lwcrv{lag}(:,2),'r','LineWidth',2);
            set(handles.length_width_axis,'Box','on','YLim',ylim);
            xlim=get(handles.length_width_axis,'XLim');
            if get(handles.checkbox1,'Value')
                line(xlim,[handles.CrR handles.CrR],'Color','g','LineWidth',2);
                line([lbase lbase],[0 handles.CrR],'Color',[.4 .4 .4]);
            end
            xlabel('length from the tip');
            ylabel('width');
            
            handles.CurrGeom=lwcrv{lag};
            handles.CurrBase=[lbase,handles.CrR];
        end    
    end
            
    handles.OLD=lag;
    guidata(hObject, handles);
   
function switchFilt(hObject, src, event) 

    handles = guidata(hObject);
    handles.CrR = round(event.NewValue);
    
    %old = handles.OLD; 
    %if handles.CrR~=old 
    
        fr = round(get(handles.movieSlider,'Value')); 
        lag = round(get(handles.lagSlider,'Value'));
        
        set(handles.CrR_text, 'String', ['CrR = ' num2str(handles.CrR,'%.2f')]);

        Tip=handles.Tip{fr};
        rpath=handles.RedPaths{fr};
        lwcrv=handles.LWcurves{fr};     

        [xbase,ybase,lbase,ibase]=FindBase(rpath{lag},lwcrv{lag},handles.CrR,handles.CutOff);

        axes(handles.boundary_axis);
        xlim=get(handles.boundary_axis,'XLim');
        ylim=get(handles.boundary_axis,'YLim');
        cla;
        hold on;
        bound=handles.BOUNDARIES{fr};
        by=bound(:,1); bx=bound(:,2); 
        fill(bx,by,[.9 .9 .9],'EdgeColor','k');
        plot(rpath{lag}(:,1),rpath{lag}(:,2),'r','LineWidth',2);
        plot(bx(Tip),by(Tip),'bo','MarkerFaceColor','b','MarkerSize',2);
        if get(handles.checkbox1,'Value')
            plot([rpath{lag}(1:ibase,1); xbase],[rpath{lag}(1:ibase,2); ybase],'g','LineWidth',2);
            plot(xbase,ybase,'go','MarkerFaceColor','r','MarkerSize',2);
        end
        axis image;
        axis ij;
        set(handles.boundary_axis,'XLim',xlim);
        set(handles.boundary_axis,'YLim',ylim);
        axis off;

        axes(handles.length_width_axis);
        xlim=get(handles.length_width_axis,'XLim');
        ylim=get(handles.length_width_axis,'YLim');
        cla reset; hold on;
        plot(lwcrv{lag}(:,1),lwcrv{lag}(:,2),'r','LineWidth',2);
        set(handles.length_width_axis,'Box','on','XLim',xlim,'YLim',ylim);    
        if get(handles.checkbox1,'Value')
            line(xlim,[handles.CrR handles.CrR],'Color','g','LineWidth',2);
            line([lbase lbase],[0 handles.CrR],'Color',[.4 .4 .4]);
        end
        xlabel('length from the tip');
        ylabel('width');

        handles.CurrGeom=lwcrv{lag};
        handles.CurrBase=[lbase,handles.CrR];
    %end
    %handles.OLD=handles.CrR;

    guidata(hObject, handles);
%--------------------------------------------------------------------------       



% --- Additional functions ------------------------------------------------
function import(hObject,imp)
    handles = guidata(hObject);
    
    handles.Lag = 1;
    handles.IMAGE_NUMBER = 1;
    handles.NeedSpeed = 0;
    try
        delete(handles.LisMov);
    end
    try
        delete(handles.LisLag);
    end
    try
        delete(handles.LisFil);
    end
    
    axes(handles.length_width_axis);
    %set(handles.length_width_axis,'Position',[550 295 500 430]);
    colormap(gray); 
    cla; axis fill;    
    set(handles.length_width_axis,'Box','on','XTick',[],'YTick',[],'Visible','on');
    
    axes(handles.boundary_axis);
    %set(handles.boundary_axis,'Position',[20 295 500 430]);
    cla; axis fill;    
    set(handles.boundary_axis,'Box','on','XTick',[],'YTick',[],'Visible','on');
    
    axes(handles.profile_axis);
    %set(handles.profile_axis,'Position',[50 35 470 180]);
    cla; axis fill;    
    set(handles.profile_axis,'Box','on','XTick',[],'YTick',[],'Visible','on');
    
    axes(handles.progress1);    
    cla reset;     
    set(handles.progress1,'Box','on','XTick',[],'YTick',[],'Visible','off');
    
    axes(handles.progress2);    
    cla reset;     
    set(handles.progress2,'Box','on','XTick',[],'YTick',[],'Visible','off');
    
    set(handles.Tip_text,'String','  ');
    set(handles.frame_text,'String','Frame','ForegroundColor',[.5 .5 .5]);
    set(handles.CrR_text,'String','  ');
    set(handles.edit_smoo,'String',' ','Enable','off');
    set(handles.results_panel,'ForegroundColor',[.5 .5 .5]);
    set(handles.smoo_text,'ForegroundColor',[.5 .5 .5]);
    set(handles.label_text,'ForegroundColor',[.5 .5 .5]);
    set(handles.bound_profile,'String','  ');      
    
    set(handles.lev_1_1,'String','  '); set(handles.lev_1_2,'String','  '); set(handles.lev_1_3,'String','  '); set(handles.lev_1_4,'String','  '); set(handles.lev_1_5,'String','  ');
    set(handles.lev_2_1,'String','  '); set(handles.lev_2_2,'String','  '); set(handles.lev_2_3,'String','  '); set(handles.lev_2_4,'String','  '); set(handles.lev_2_5,'String','  ');
    set(handles.lev_3_1,'String','  '); set(handles.lev_3_2,'String','  '); set(handles.lev_3_3,'String','  '); set(handles.lev_3_4,'String','  '); set(handles.lev_3_5,'String','  ');
    
    set(handles.movieSlider, 'Enable', 'off');
    set(handles.lagSlider, 'Enable', 'off');
    set(handles.CrR_Slider, 'Enable', 'off');
        
    set(handles.graph_frame_button, 'Enable', 'off');
    set(handles.graph_all_button, 'Enable', 'off');
    set(handles.profile_button, 'Enable', 'off');
    set(handles.batch_button, 'Enable', 'off');
    set(handles.list_button, 'Enable', 'off');
    
    set(handles.uitable1,'Data',cell(1,2),'Enable','off');
    
    set(handles.checkbox1, 'Value', 0);
    set(handles.checkbox1, 'Enable', 'off');
    
    set(handles.save_res, 'Enable', 'off');
    
    if imp==1
        [filename,pathname,filterindex] = uigetfile('*.mat'); 
    elseif imp==2
        [filename,pathname,filterindex] = uigetfile('*.tif'); 
    end
    
    if filename~=0
        
        jFigPeer = get(handle(gcf),'JavaFrame');
        try
            jWindow = jFigPeer.fFigureClient.getWindow;
        catch
            jWindow = jFigPeer.fHG1Client.getWindow;
        end
        set(handle(jWindow),'Enabled',false);

        try

            okey=0;

            if imp==1

                bnd=load([pathname filename],'boundaries');

                if ~isempty(bnd)
                    okey=1;


                    b=bnd.boundaries;
                    handles.NUM_ELEMENTS=length(b);

                    xmi= zeros(1,handles.NUM_ELEMENTS);
                    xma=xmi; ymi=xmi; yma=xmi; 
                    axes(handles.progress1);
                    cla;
                    for k = 1:handles.NUM_ELEMENTS               
                       bx=[0 k k 0 0]/handles.NUM_ELEMENTS;
                       by=[0 0 1 1 0];
                       fill(bx,by,[.7 .7 .7]);
                       set(handles.progress1,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[]);
                       pause(0.001);
                       %tmp=round(b{k}); 
                       [tmp,~]=regbound(b{k});
                       b{k}=tmp;                                  
                       xmi(k)=min(tmp(:,2)); xma(k)=max(tmp(:,2));
                       ymi(k)=min(tmp(:,1)); yma(k)=max(tmp(:,1));                   
                    end
                    cla(handles.progress1);
                    set(handles.progress1,'Visible','off');

                    xlim=[min(xmi)-10 max(xma)+10];
                    ylim=[min(ymi)-10 max(yma)+10];               

                end

            elseif imp==2
                info = imfinfo([pathname filename]);
                handles.NUM_ELEMENTS = numel(info);
                if handles.NUM_ELEMENTS>0
                    okey=1;

                    b = cell(1,handles.NUM_ELEMENTS);
                    xmi= zeros(1,handles.NUM_ELEMENTS);
                    xma=xmi; ymi=xmi; yma=xmi; 

                    axes(handles.progress1);
                    cla;
                    for k = 1:handles.NUM_ELEMENTS               
                        bx=[0 k k 0 0]/handles.NUM_ELEMENTS;
                        by=[0 0 1 1 0];
                        fill(bx,by,[.7 .7 .7]);
                        set(handles.progress1,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[]);
                        pause(0.001);                    
                        image = imread([pathname filename], k);
                        if length(size(image))==3
                            image=image(:,:,1);
                        end
                        image = cast(image, 'double');
                        image(image>0) = 1;
                        if max(image(:))==min(image(:))
                            image=ones(size(image));
                        end

                        [B,L]=bwboundaries(image,'noholes');
                        A=regionprops(logical(image),'Area');
                        Area=[A.Area];
                        [~,ind]=max(Area);
                        if ~isempty(ind)
                            %tmp=B{ind}; 
                            [tmp,~]=regbound(B{ind});
                            L(L~=ind)=0; L(L==ind)=1;
                            smx=sum(L,1); smy=sum(L,2);
                            xmi(k)=find(smx>0,1,'first');
                            xma(k)=find(smx>0,1,'last');
                            ymi(k)=find(smy>0,1,'first');
                            yma(k)=find(smy>0,1,'last');
                            b{k}=tmp;                                                                                              
                        else
                            xmi(k)=Inf;
                            xma(k)=-Inf;
                            ymi(k)=Inf;
                            yma(k)=-Inf;
                            %lng(k)=1;
                        end

                    end
                    cla(handles.progress1);
                    set(handles.progress1,'Visible','off');

                    xlim=[min(xmi)-10 max(xma)+10];
                    ylim=[min(ymi)-10 max(yma)+10];

                end            
            end

            if okey


                handles.GN = handles.DefGN;
                handles.CrR = handles.DefCrR;

                handles.BOUNDARIES = b;
                handles.X = xlim;
                handles.Y = ylim;

                handles.REM1=cell(1,handles.NUM_ELEMENTS);
                handles.REM2=cell(1,handles.NUM_ELEMENTS);
                handles.Nconvex=zeros(1,handles.NUM_ELEMENTS);
                handles.Nconcave=zeros(1,handles.NUM_ELEMENTS);
                handles.Nedge=zeros(1,handles.NUM_ELEMENTS);  
                
                handles.gbound=cell(1,handles.NUM_ELEMENTS);  

                handles.INDX=cell(1,handles.NUM_ELEMENTS);
                handles.Vb=cell(1,handles.NUM_ELEMENTS);
                handles.Coor=cell(1,handles.NUM_ELEMENTS);
                handles.DIST=cell(1,handles.NUM_ELEMENTS);
                handles.Rad=cell(1,handles.NUM_ELEMENTS);
                handles.Pr=cell(1,handles.NUM_ELEMENTS);
                handles.Py=cell(1,handles.NUM_ELEMENTS);
                handles.Tip=cell(1,handles.NUM_ELEMENTS);

                handles.PathInd=cell(1,handles.NUM_ELEMENTS);
                handles.PathLng=cell(1,handles.NUM_ELEMENTS);
                                
                handles.DoneGR=zeros(1,handles.NUM_ELEMENTS);
                handles.DonePR=zeros(1,handles.NUM_ELEMENTS);

                axes(handles.boundary_axis);
                hold on;
                bound=b{1};
                by=bound(:,1); bx=bound(:,2);
                fill(bx,by,[.9 .9 .9],'EdgeColor','k');  
                axis image; 
                axis ij;
                set(handles.boundary_axis,'XLim',xlim,'YLim',ylim);
                axis off;

                if handles.NUM_ELEMENTS>1
                    set(handles.movieSlider, 'Min', 1, 'Max', handles.NUM_ELEMENTS,'Value',1);
                    set(handles.movieSlider, 'SliderStep', [1/(handles.NUM_ELEMENTS-1) 5/(handles.NUM_ELEMENTS-1)]);
                    set(handles.movieSlider, 'Enable', 'on');
                else
                    set(handles.movieSlider, 'Value',1);
                end
                %set(handles.lagSlider,'Min', 1, 'Max', handles.NUM_ELEMENTS-1,'Value',1);

                set(handles.frame_text,'ForegroundColor','k');
                set(handles.frame_text,'String','Frame # 1');

                set(handles.results_panel,'ForegroundColor','k');
                set(handles.label_text,'ForegroundColor','k');                        

                set(handles.graph_frame_button, 'Enable', 'on');
                set(handles.graph_all_button, 'Enable', 'on');  
                set(handles.list_button, 'Enable', 'on');

                

                [namestr,handles.nm1,handles.nm2,handles.nm3,handles.ftype]=batchfolders(handles,[pathname filename]);
                sdata=cell(1,2);
                sdata{1,1}=namestr; sdata{1,2}='    imported';
                set(handles.uitable1,'Data',sdata,'Enable','on');

                guidata(hObject, handles);
            end
        catch
            errordlg('Import failed. Please check the file format.');    
        end
        set(handle(jWindow),'Enabled',true);
    end
    

    handles.LisMov = addlistener(handles.movieSlider, 'Value', 'PostSet', @(src, event)switchFrame(hObject, src, event));
    handles.LisLag = addlistener(handles.lagSlider, 'Value', 'PostSet', @(src, event)switchLag(hObject, src, event));
    handles.LisFil = addlistener(handles.CrR_Slider, 'Value', 'PostSet', @(src, event)switchFilt(hObject, src, event));

    guidata(hObject,handles);
     
function handles=replot(handles)

    fr = round(get(handles.movieSlider,'Value')); 
    lag = round(get(handles.lagSlider,'Value'));  
    %handles.CrR = round(get(handles.lagSlider,'Value'));  
        
    axes(handles.boundary_axis);
    xlim=get(handles.boundary_axis,'XLim');
    ylim=get(handles.boundary_axis,'YLim');
    cla;
    hold on;
            
    bound=handles.BOUNDARIES{fr};
    by=bound(:,1); bx=bound(:,2); 
    fill(bx,by,[.9 .9 .9],'EdgeColor','k');
    axis image;
    axis ij;
    set(handles.boundary_axis,'XLim',xlim);
    set(handles.boundary_axis,'YLim',ylim);
    axis off;
    
    if handles.DoneGR(fr)
        REM1=handles.REM1{fr};
        REM2=handles.REM2{fr};
        
        axes(handles.boundary_axis);
        hold on;
        plot([REM2(:,1)'; REM2(:,3)'],[REM2(:,2)';REM2(:,4)'],'Color',[.8 .8 .8]);
        plot([REM1(:,1)'; REM1(:,3)'],[REM1(:,2)';REM1(:,4)'],'Color',[.2 .2 .2]);        
        
        if handles.DonePR(fr)             
                       
            Yb=by; Xb=bx;    
            INDX=handles.INDX{fr};
            Vb=handles.Vb{fr};
            Coor=handles.Coor{fr};
            DIST=handles.DIST{fr};            
            Tip=handles.Tip{fr};
            Pr=handles.Pr{fr};
            Py=handles.Py{fr};               
            Rad=handles.Rad{fr};     
            
            PathInd=handles.PathInd{fr};
            PathLng=handles.PathLng{fr};
            
            rpath=handles.RedPaths{fr};
            %gpath=handles.GrnPaths{fr};  %might not need
            lwcrv=handles.LWcurves{fr}; 
            
            %if handles.CrR>2*max(Rad);
            %    handles.CrR=2*max(Rad);
            %end
            %if lag>length(Tip)
            %    lag=length(Tip);
            %end
            
            [xbase,ybase,lbase,ibase]=FindBase(rpath{lag},lwcrv{lag},handles.CrR,handles.CutOff);
            
            if get(handles.checkbox1,'Value')
                Lt=length(Tip);
                tuck=zeros(1,Lt);
                for i=1:Lt
                    [~,filo,~]=FindPathCond(INDX,Coor,DIST,Rad,Vb(Tip(i)),handles.CrR,handles.CutOff);
                    if filo(end,1)>0
                        tuck(i)=filo(end,1);
                    else
                        tuck(i)=filo(end-1,1);
                    end
                end

                nXb=zeros(size(Xb));
                nYb=nXb; nLb=0;
                for i=1:length(Xb)
                    tpath=PathInd(1:PathLng(i),i);
                    for j=1:Lt
                        tmp=find(tpath==tuck(j),1);
                        if ~isempty(tmp)
                            break;
                        end
                    end
                    if isempty(tmp)
                        nLb=nLb+1;
                        nXb(nLb)=Xb(i);
                        nYb(nLb)=Yb(i);
                    end
                end
                cXb=nXb(1:nLb);
                cYb=nYb(1:nLb);

                if ~isempty(cXb)
                    handles.gbound{fr}=[[cYb; cYb(1)],[cXb; cXb(1)]];
                    %axes(handles.boundary_axis); 
                    plot([cXb; cXb(1)]',[cYb; cYb(1)]','Color','g','LineWidth',2);
                else
                    handles.gbound{fr}=[];
                end
            end
            plot(Xb(Tip),Yb(Tip),'bo','MarkerFaceColor','b','MarkerSize',2);            
            plot(rpath{lag}(:,1),rpath{lag}(:,2),'r','LineWidth',2);
            
                       
            axes(handles.profile_axis);
            cla reset;
            hold on;
            plot(Pr,'c');
            plot(Py,'k','LineWidth',1);
            plot(Tip,Py(Tip),'bo','MarkerFaceColor','b','MarkerSize',4);
            set(handles.profile_axis,'Box','on','XLim',[-2 length(Xb)+2]);            
            ylim=get(handles.profile_axis,'YLim');            
            line([Tip(lag) Tip(lag)],ylim,'Color','r','LineWidth',2);
            
            axes(handles.length_width_axis);
            cla reset; hold on;
            plot(lwcrv{lag}(:,1),lwcrv{lag}(:,2),'r','LineWidth',2);
            set(handles.length_width_axis,'Box','on','YLim',[0 max(Rad)]);
            xlim=get(handles.length_width_axis,'XLim');
            if get(handles.checkbox1,'Value')
                line(xlim,[handles.CrR handles.CrR],'Color','g','LineWidth',2);
                line([lbase lbase],[0 handles.CrR],'Color',[.4 .4 .4]);
            end
            xlabel('length from the tip');
            ylabel('width');            
            
            handles.CurrGeom=lwcrv{lag};
            handles.CurrBase=[lbase,handles.CrR];
        end       
    end    
    
function [namestr,nm1,nm2,nm3,ftype]=batchfolders(handles,str)

    subdir=cell(20,1);
    n=zeros(20,1);
    k=1; 
    for i=1:length(str)
        if double(str(i))~=92 && double(str(i))~=47 
            subdir{k}=[subdir{k} str(i)];
            n(k)=n(k)+1;
        else
            k=k+1;
        end
    end

    subdir=subdir(1:k);
    ftype=subdir{k}((n(k)-2):n(k));
    subdir{k}=subdir{k}(1:(n(k)-4));
    digdir=cell(k,1);
    remdir=cell(k,1);
    rr=zeros(k,1);
    dg=zeros(k,1);
    n=n(1:k);
    n(k)=n(k)-4;
    for i=1:k
        str=subdir{i};
        for j=n(i):-1:1
           if double(str(j))>47 && double(str(j))<58
               remdir{i}=[str(j) remdir{i}];
               continue;
           else
               digdir{i}=str(1:j);
               break;
           end
        end 
        if ~isempty(remdir{i});
            rr(i)=str2double(remdir{i});
            dg(i)=length(remdir{i});
        end
    end
    
    if k>3
        
        nm1=subdir{1};
        for q=2:(k-3)
            nm1=[nm1 '/' subdir{q}];
        end
        nm1=[nm1 '/' digdir{k-2}];
        nm2=digdir{k-1}; nm3=digdir{k};
        namestr=['../' subdir{k-2} '/' subdir{k-1} '/' subdir{k} '.' ftype];
    elseif k>2
        nm1=[];
        nm2=subdir{1};
        for q=2:(k-2)
            nm2=[nm2 '/' subdir{q}];
        end
        nm2=[nm2 '/' digdir{k-1}];
        nm3=digdir{k};
        namestr=['../' subdir{k-1} '/' subdir{k} '.' ftype];
    elseif k>1
        nm1=[]; nm2=[]; 
        nm3=subdir{1};
        for q=2:(k-1)
            nm3=[nm3 '/' subdir{q}];
        end        
        nm3=[nm3 '/' digdir{k}];
        namestr=['../' subdir{k} '.' ftype];
    end
    
        
    if k>1   
        set(handles.lev_1_1,'String',digdir{k}); 
        if rr(k)>0
            set(handles.lev_1_2,'String',rr(k));
            set(handles.lev_1_3,'String',rr(k)); 
        end
        set(handles.lev_1_4,'String','1');
        set(handles.lev_1_5,'String',num2str(dg(k)));
    end
    
    if k>2
        set(handles.lev_2_1,'String',digdir{k-1}); 
        if rr(k-1)>0
            set(handles.lev_2_2,'String',rr(k-1));
            set(handles.lev_2_3,'String',rr(k-1)); 
        end
        set(handles.lev_2_4,'String','1');
        set(handles.lev_2_5,'String',num2str(dg(k-1)));

    end
    if k>3
        set(handles.lev_3_1,'String',digdir{k-2}); 
        if rr(k-2)>0
            set(handles.lev_3_2,'String',rr(k-2));
            set(handles.lev_3_3,'String',rr(k-2)); 
        end
        set(handles.lev_3_4,'String','1');
        set(handles.lev_3_5,'String',num2str(dg(k-2)));

    end    
%--------------------------------------------------------------------------    
  


% --- Sliders -------------------------------------------------------------
function movieSlider_Callback(hObject, eventdata, handles)

    set(handles.movieSlider,'Enable','off');
    set(handles.lagSlider,'Enable','off');
    set(handles.CrR_Slider,'Enable','off');
    
    handles = guidata(hObject); 
    fr = round(get(handles.movieSlider,'Value')); 
    lag = round(get(handles.lagSlider,'Value'));
    %set(handles.movieSlider,'Value',fr);
    set(handles.frame_text,'String',['Frame # ' num2str(fr)]);   

    if handles.DoneGR(fr)
       handles=replot(handles); 
    end
    
    if handles.DonePR(fr) 
        
        handles.OLD = lag;
        
        Rad=handles.Rad{fr};
        Tip=handles.Tip{fr};
        
        if handles.CrR>floor(max(Rad))
            handles.CrR=floor(max(Rad));
        end
        if lag>length(Tip)
            lag=1;%length(Tip);
        end
        
        try
            delete(handles.LisFil);
        end
        set(handles.CrR_Slider,'Min',1,'Max',floor(max(Rad)),'Value',handles.CrR);
        set(handles.CrR_Slider,'SliderStep',[1/(floor(max(Rad))-1) 5/(floor(max(Rad))-1)]);                
        set(handles.CrR_text, 'String', ['CrR = ' num2str(handles.CrR,'%.2f')]);            
        handles.LisFil = addlistener(handles.CrR_Slider, 'Value', 'PostSet', @(src, event)switchFilt(hObject, src, event));

        try
            delete(handles.LisLag);
        end
        set(handles.lagSlider,'Min',1,'Max',length(Tip),'Value',lag);
        set(handles.lagSlider,'SliderStep',[1/(length(Tip)-1) 5/(length(Tip)-1)]);
        set(handles.Tip_text, 'String', ['Tip # ' num2str(lag)]);
        handles.LisLag = addlistener(handles.lagSlider, 'Value', 'PostSet', @(src, event)switchLag(hObject, src, event));        
        
        handles=replot(handles);
    end
    
    if handles.NUM_ELEMENTS>1
        set(handles.movieSlider,'Enable','on');
    end
    set(handles.lagSlider,'Enable','on');
    if get(handles.checkbox1,'Value')
        set(handles.CrR_Slider,'Enable','on');
    end
    
    guidata(hObject, handles);

function lagSlider_Callback(hObject, eventdata, handles)

    set(handles.movieSlider,'Enable','off');
    set(handles.lagSlider,'Enable','off');
    set(handles.CrR_Slider,'Enable','off');
    
    handles = guidata(hObject);     
    lag = round(get(handles.lagSlider,'Value'));
    %set(handles.lagSlider,'Value',lag);
    set(handles.Tip_text, 'String', ['Tip # ' num2str(lag)]);
    
    handles=replot(handles);
       
    %set(handles.addtoqueue,'Enable','on');
    if handles.NUM_ELEMENTS>1
        set(handles.movieSlider,'Enable','on');
    end
    set(handles.lagSlider,'Enable','on');
    if get(handles.checkbox1,'Value')
        set(handles.CrR_Slider,'Enable','on');
    end
    guidata(hObject, handles);
    
function CrR_Slider_Callback(hObject, eventdata, handles)

    set(handles.movieSlider,'Enable','off');
    set(handles.lagSlider,'Enable','off');
    set(handles.CrR_Slider,'Enable','off');
    
    handles = guidata(hObject);     
    %handles.CrR = round(get(handles.CrR_Slider,'Value'));
    set(handles.CrR_text, 'String', ['CrR = ' num2str(handles.CrR,'%.2f')]);
    
    handles=replot(handles);
       
    %set(handles.addtoqueue,'Enable','on');
    if handles.NUM_ELEMENTS>1
        set(handles.movieSlider,'Enable','on');
    end
    set(handles.lagSlider,'Enable','on');
    set(handles.CrR_Slider,'Enable','on');
    guidata(hObject, handles);
    
function movieSlider_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end

function lagSlider_CreateFcn(hObject, eventdata, handles)

    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
    
function CrR_Slider_CreateFcn(hObject, eventdata, handles)

    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
%--------------------------------------------------------------------------

    
  
% --- Buttons -------------------------------------------------------------
function graph_frame_button_Callback(hObject, eventdata, handles)

jFigPeer = get(handle(gcf),'JavaFrame');
try
    jWindow = jFigPeer.fFigureClient.getWindow;
catch
    jWindow = jFigPeer.fHG1Client.getWindow;
end
set(handle(jWindow),'Enabled',false);

try
    
    fr=handles.IMAGE_NUMBER;
    
    if sum(handles.DoneGR)==0
        axes(handles.progress2);   
        cla reset; 
        set(handles.progress2,'Box','on','XTick',[],'YTick',[],'Visible','on');
    end
    
    if handles.DoneGR(fr)
        %REM1=handles.REM1{fr};
        %REM2=handles.REM2{fr};
        %axes(handles.boundary_axis);
        %hold on;
        %plot([REM2(:,1)'; REM2(:,3)'],[REM2(:,2)';REM2(:,4)'],'Color',[.8 .8 .8]);
        %plot([REM1(:,1)'; REM1(:,3)'],[REM1(:,2)';REM1(:,4)'],'Color',[.2 .2 .2]);
        set(handle(jWindow),'Enabled',true);
        return;
    end
    
    bound=round(handles.BOUNDARIES{fr});
    by=bound(:,1); bx=bound(:,2); 
    
    Status=0; atmp=0; tm=0;
    
    while ~Status && atmp<30
        atmp=atmp+1;       
            
        Xb=bx+0.02*(randn(size(bx))-0.5);
        Yb=by+0.02*(randn(size(bx))-0.5);
    
        axes(handles.boundary_axis);
        xlim=get(handles.boundary_axis,'XLim');
        ylim=get(handles.boundary_axis,'YLim');
        cla;
        hold on;
        fill(Xb,Yb,[.9 .9 .9],'EdgeColor','k');
        axis image;
        axis ij;
        set(handles.boundary_axis,'XLim',xlim);
        set(handles.boundary_axis,'YLim',ylim);
        axis off;

        try
            tic
            [REM1,REM2,Status,Nconvex,Nconcave,Nedge]=FindGraphX(Xb',Yb',length(bx),handles);
            tm=tm+toc;
        catch
            Status=0;
        end
    end
    disp(' ');
    disp(['Frame # ' num2str(fr) '; Attempts ' num2str(atmp) '; took ' num2str(round(tm)) ' sec.']);
    
    if Status
        
        axes(handles.boundary_axis);
        hold on;
        plot([REM2(:,1)'; REM2(:,3)'],[REM2(:,2)';REM2(:,4)'],'Color',[.8 .8 .8]);
        plot([REM1(:,1)'; REM1(:,3)'],[REM1(:,2)';REM1(:,4)'],'Color',[.2 .2 .2]);
        
        handles.DoneGR(fr)=1;
        
        handles.REM1{fr}=REM1;
        handles.REM2{fr}=REM2;
        handles.Nconvex(fr)=Nconvex;
        handles.Nconcave(fr)=Nconcave;
        handles.Nedge(fr)=Nedge;  
        handles.BOUNDARIES{fr}=[Yb,Xb];
        
        axes(handles.progress2);
        cla;
        hold on;
        k=sum(handles.DoneGR);
        bx=[0 k k 0 0]/handles.NUM_ELEMENTS;
        by=[0 0 1 1 0];
        fill(bx,by,[.7 .7 .7]);
        set(handles.progress2,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[]);
        
        set(handles.profile_button, 'Enable', 'on');        
                
    end
    guidata(hObject, handles); 
    
end
set(handle(jWindow),'Enabled',true);
 
function graph_all_button_Callback(hObject, eventdata, handles)

jFigPeer = get(handle(gcf),'JavaFrame');
try
    jWindow = jFigPeer.fFigureClient.getWindow;
catch
    jWindow = jFigPeer.fHG1Client.getWindow;
end
set(handle(jWindow),'Enabled',false);

try

    if sum(handles.DoneGR)==0
        axes(handles.progress2);
        cla reset;    
        set(handles.progress2,'Box','on','XTick',[],'YTick',[],'Visible','on');
    end
    
    for fr=1:handles.NUM_ELEMENTS 
        
        set(handles.movieSlider,'Value',fr);
        handles=guidata(hObject);
        if ~handles.DoneGR(fr)
            graph_frame_button_Callback(hObject, eventdata, handles);
            handles.DoneGR(fr)=1;            
        end
        
        axes(handles.progress2);
        cla;
        hold on;
        k=sum(handles.DoneGR);
        bx=[0 k k 0 0]/handles.NUM_ELEMENTS;
        by=[0 0 1 1 0];
        fill(bx,by,[.7 .7 .7]);
        set(handles.progress2,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[]);
        pause(0.001);                        
    end    
    handles=guidata(hObject);
    guidata(hObject, handles);  
 
end
set(handle(jWindow),'Enabled',true);    
        
function profile_button_Callback(hObject, eventdata, handles)

jFigPeer = get(handle(gcf),'JavaFrame');
try
    jWindow = jFigPeer.fFigureClient.getWindow;
catch
    jWindow = jFigPeer.fHG1Client.getWindow;
end
set(handle(jWindow),'Enabled',false);

try

    oldV=handles.IMAGE_NUMBER;
    set(handles.profile_button, 'Enable', 'off');

    for fr=1:handles.NUM_ELEMENTS
        if ~handles.DonePR(fr) && handles.DoneGR(fr)        
            
            REM=[handles.REM1{fr}; handles.REM2{fr}];
            Smoo=handles.GN;

            bound=handles.BOUNDARIES{fr};
            Yb=bound(:,1); Xb=bound(:,2);

            Nconvex=handles.Nconvex(fr);
            Nconcave=handles.Nconcave(fr);
            Nedge=handles.Nedge(fr);  

            [INDX,Vb,Coor,DIST1,DIST2,Rad,Pr,Py,Tip]=FindProfile(REM,Xb,Yb,Nconvex,Nconcave,Nedge,Smoo);            
            [PathInd,PathLng]=PathIndex(INDX,Vb);
            DIST=(DIST1+DIST2)/2;                        
            
            %---------------------------------
            Lt=length(Tip);
            tuck=zeros(1,Lt);
            rpath=cell(1,Lt);
            gpath=cell(1,Lt);
            lwcrv=cell(1,Lt);
            
            for i=1:Lt                
                [path,filo,~]=FindPathCond(INDX,Coor,DIST,Rad,Vb(Tip(i)),handles.CrR,handles.CutOff);
                
                tpath=PathInd(1:PathLng(Tip(i)),Tip(i));
                if filo(end,1)>0
                    tuck(i)=filo(end,1);
                else
                    tuck(i)=filo(end-1,1);
                end
                tmp=find(tpath==tuck(i),1);                
                
               rpath{i}=[[path(:,2); Coor(INDX(1,2),1)],[path(:,3); Coor(INDX(1,2),2)]];
               if ~isempty(tmp)
                   gpath{i}=[[Coor(tpath(1:tmp),1); filo(end,2)],[Coor(tpath(1:tmp),2); filo(end,3)]];                
               end
               lwcrv{i}=[[max(path(:,4))-path(:,4); max(path(:,4))],[path(:,5); max(Rad)]];
               
            end
            
            nXb=zeros(size(Xb));
            nYb=nXb; nLb=0;
            for i=1:length(Xb)
                tpath=PathInd(1:PathLng(i),i);
                for j=1:Lt
                    tmp=find(tpath==tuck(j),1);
                    if ~isempty(tmp)
                        break;
                    end
                end
                if isempty(tmp)
                    nLb=nLb+1;
                    nXb(nLb)=Xb(i);
                    nYb(nLb)=Yb(i);
                end
            end
            cXb=nXb(1:nLb);
            cYb=nYb(1:nLb);

            if ~isempty(cXb)
                handles.gbound{fr}=[[cYb; cYb(1)],[cXb; cXb(1)]];
                %axes(handles.boundary_axis); 
                plot([cXb; cXb(1)]',[cYb; cYb(1)]','Color','g','LineWidth',2);
            else
                handles.gbound{fr}=[];
            end
            %----------------------------------
            
            handles.DonePR(fr)=1;
                        
            axes(handles.progress2);
            cla;
            hold on;
            k=sum(handles.DonePR);
            bx=[0 k k 0 0]/handles.NUM_ELEMENTS;
            by=[0 0 1 1 0];
            fill(bx,by,[.7 .7 .7]);
            set(handles.progress2,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[]);
            pause(0.001);     

            handles.INDX{fr}=INDX;
            handles.Vb{fr}=Vb;
            handles.Coor{fr}=Coor;
            handles.DIST{fr}=DIST;
            handles.Pr{fr}=Pr;
            handles.Py{fr}=Py;
            handles.Rad{fr}=Rad;
            handles.Tip{fr}=Tip;
            
            handles.PathInd{fr}=PathInd;
            handles.PathLng{fr}=PathLng;                
            
            handles.RedPaths{fr}=rpath;
            %handles.GrnPaths{fr}=gpath;
            handles.LWcurves{fr}=lwcrv;                                
            %handles.TipCoord{fr}=[Xb(Tip),Yb(Tip)];

            guidata(hObject, handles); 
        end
    end
    
    set(handles.bound_profile,'String','Boundary Profile');
    
    set(handles.checkbox1, 'Enable', 'on');
    set(handles.lagSlider, 'Enable', 'on');
    set(handles.edit_smoo,'String',num2str(handles.GN),'Enable','on');
    set(handles.smoo_text,'ForegroundColor','k');
    %set(handles.mupx_text,'ForegroundColor','k');
    %set(handles.scale_text,'ForegroundColor','k');
    
    set(handles.save_res, 'Enable', 'on');
    %set(handles.addtoqueue, 'Enable', 'on');
    %set(handles.del_sel, 'Enable', 'on');
    set(handles.uitable1,'Enable','on');
    %set(handles.scale_edit, 'Enable', 'on');
    
    movieSlider_Callback(hObject, eventdata, handles);  
    
end
set(handle(jWindow),'Enabled',true);    

function list_button_Callback(hObject, eventdata, handles)

jFigPeer = get(handle(gcf),'JavaFrame');
try
    jWindow = jFigPeer.fFigureClient.getWindow;
catch
    jWindow = jFigPeer.fHG1Client.getWindow;
end
set(handle(jWindow),'Enabled',false);

try

    handles.fnames=cell(1,1000);
    fls=0;

    if isempty(handles.nm1) && isempty(handles.nm2)

        st3=str2num(get(handles.lev_1_2,'String'));
        en3=str2num(get(handles.lev_1_3,'String'));
        in3=str2num(get(handles.lev_1_4,'String'));
        dg3=str2num(get(handles.lev_1_5,'String'));

        leng=length(get(handles.lev_1_1,'String'))+dg3+4;

        if ~isempty(st3) && ~isempty(en3) && dg3>0
            frm3=['%0' num2str(dg3) 'd'];
            for i=st3:in3:en3
                fls=fls+1;
                handles.fnames{fls}=[handles.nm3 num2str(i,frm3) '.' handles.ftype];            
            end
        else
            fls=fls+1;
            handles.fnames{fls}=[handles.nm3 '.' handles.ftype];
        end
    elseif isempty(handles.nm1) && ~isempty(handles.nm2)    

        st3=str2num(get(handles.lev_1_2,'String'));
        en3=str2num(get(handles.lev_1_3,'String'));
        in3=str2num(get(handles.lev_1_4,'String'));
        dg3=str2num(get(handles.lev_1_5,'String'));

        st2=str2num(get(handles.lev_2_2,'String'));
        en2=str2num(get(handles.lev_2_3,'String'));
        in2=str2num(get(handles.lev_2_4,'String'));
        dg2=str2num(get(handles.lev_2_5,'String'));

        leng=length(get(handles.lev_1_1,'String'))+length(get(handles.lev_2_1,'String'))+dg2+dg3+5;    

        if ~isempty(st2) && ~isempty(en2) && dg2>0
            frm2=['%0' num2str(dg2) 'd'];
            for i=st2:in2:en2            
                if ~isempty(st3) && ~isempty(en3) && dg3>0
                    frm3=['%0' num2str(dg3) 'd'];
                    for j=st3:in3:en3
                        fls=fls+1;
                        handles.fnames{fls}=[handles.nm2 num2str(i,frm2) '/' handles.nm3 num2str(j,frm3) '.' handles.ftype];
                    end
                else
                    fls=fls+1;
                    handles.fnames{fls}=[handles.nm2 num2str(i,frm2) '/' handles.nm3 '.' handles.ftype];
                end
            end
        else
            if ~isempty(st3) && ~isempty(en3) && dg3>0
                frm3=['%0' num2str(dg3) 'd'];
                for j=st3:in3:en3
                    fls=fls+1;
                    handles.fnames{fls}=[handles.nm2 '/' handles.nm3 num2str(j,frm3) '.' handles.ftype];
                end
            else
                fls=fls+1;
                handles.fnames{fls}=[handles.nm2 '/' handles.nm3 '.' handles.ftype];
            end
        end


    elseif ~isempty(handles.nm1) && ~isempty(handles.nm2)    

        st3=str2num(get(handles.lev_1_2,'String'));
        en3=str2num(get(handles.lev_1_3,'String'));
        in3=str2num(get(handles.lev_1_4,'String'));
        dg3=str2num(get(handles.lev_1_5,'String'));

        st2=str2num(get(handles.lev_2_2,'String'));
        en2=str2num(get(handles.lev_2_3,'String'));
        in2=str2num(get(handles.lev_2_4,'String'));
        dg2=str2num(get(handles.lev_2_5,'String'));

        st1=str2num(get(handles.lev_3_2,'String'));
        en1=str2num(get(handles.lev_3_3,'String'));
        in1=str2num(get(handles.lev_3_4,'String'));
        dg1=str2num(get(handles.lev_3_5,'String'));

        leng=length(get(handles.lev_1_1,'String'))+length(get(handles.lev_2_1,'String'))...
            +length(get(handles.lev_3_1,'String'))+dg1+dg2+dg3+6;

        if ~isempty(st1) && ~isempty(en1) && dg1>0
            frm1=['%0' num2str(dg1) 'd'];
            for i=st1:in1:en1            
                if ~isempty(st2) && ~isempty(en2) && dg2>0
                    frm2=['%0' num2str(dg2) 'd'];
                    for j=st2:in2:en2            
                        if ~isempty(st3) && ~isempty(en3) && dg3>0
                            frm3=['%0' num2str(dg3) 'd'];
                            for k=st3:in3:en3
                                fls=fls+1;
                                handles.fnames{fls}=[handles.nm1 num2str(i,frm1) '/' handles.nm2 num2str(j,frm2) '/' handles.nm3 num2str(k,frm3) '.' handles.ftype];
                            end
                        else
                            fls=fls+1;
                            handles.fnames{fls}=[handles.nm1 num2str(i,frm1) '/' handles.nm2 num2str(j,frm2) '/' handles.nm3 '.' handles.ftype];
                        end
                    end
                else
                    if ~isempty(st3) && ~isempty(en3) && dg3>0
                        frm3=['%0' num2str(dg3) 'd'];
                        for k=st3:in3:en3
                            fls=fls+1;
                            handles.fnames{fls}=[handles.nm1 num2str(i,frm1) '/' handles.nm2 '/' handles.nm3 num2str(k,frm3) '.' handles.ftype];
                        end
                    else
                        fls=fls+1;
                        handles.fnames{fls}=[handles.nm1 num2str(i,frm1) '/' handles.nm2 '/' handles.nm3 '.' handles.ftype];
                    end
                end   
            end
        else
            if ~isempty(st2) && ~isempty(en2) && dg2>0
                frm2=['%0' num2str(dg2) 'd'];
                for j=st2:in2:en2            
                    if ~isempty(st3) && ~isempty(en3) && dg3>0
                        frm3=['%0' num2str(dg3) 'd'];
                        for k=st3:in3:en3
                            fls=fls+1;
                            handles.fnames{fls}=[handles.nm1 '/' handles.nm2 num2str(j,frm2) '/' handles.nm3 num2str(k,frm3) '.' handles.ftype];
                        end
                    else
                        fls=fls+1;
                        handles.fnames{fls}=[handles.nm1 '/' handles.nm2 num2str(j,frm2) '/' handles.nm3 '.' handles.ftype];
                    end
                end
            else
                if ~isempty(st3) && ~isempty(en3) && dg3>0
                    frm3=['%0' num2str(dg3) 'd'];
                    for k=st3:in3:en3
                        fls=fls+1;
                        handles.fnames{fls}=[handles.nm1 '/' handles.nm2 '/' handles.nm3 num2str(k,frm3) '.' handles.ftype];
                    end
                else
                    fls=fls+1;
                    handles.fnames{fls}=[handles.nm1 '/' handles.nm2 '/' handles.nm3 '.' handles.ftype];
                end
            end   
        end
    end

    lenh=length(get(handles.lev_1_1,'String'))+dg3+4;
    handles.fnames=handles.fnames(1:fls);

    fnd=0;

    if fls>0
        efile=cell(1,fls);
        sdata=cell(fls,2);
        %set(handles.uitable1,'Data',sdata,'Enable','on');
        for ff=1:fls


            if exist(handles.fnames{ff},'file')==2
                fnd=fnd+1;
                efile{fnd}=handles.fnames{ff};
                sdata{fnd,1}=['..' handles.fnames{ff}((end-leng):end)];
                sdata{fnd,2}=' file is ready';
            %else
            %    sdata{ff,2}='file is not found';        

            end   

        end
        sdata=sdata(1:fnd,:);
        efile=efile(1:fnd);
    else
        sdata=cell(fls,2);
    end
    set(handles.uitable1,'Data',sdata);
    
    handles.fnd=fnd;
    handles.sdata=sdata;
    handles.efile=efile;
    handles.dg3=dg3;
    handles.lenh=lenh;
    
    set(handles.batch_button, 'Enable', 'on');
    
    guidata(hObject, handles);
    
end
set(handle(jWindow),'Enabled',true);

function batch_button_Callback(hObject, eventdata, handles)

jFigPeer = get(handle(gcf),'JavaFrame');
try
    jWindow = jFigPeer.fFigureClient.getWindow;
catch
    jWindow = jFigPeer.fHG1Client.getWindow;
end
set(handle(jWindow),'Enabled',false);

try
    
    cla(handles.progress1);
    set(handles.progress1,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[],'Visible','on');
    cla(handles.progress2);
    set(handles.progress2,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[],'Visible','on');
    
    
    fnd=handles.fnd;
    sdata=handles.sdata;
    efile=handles.efile;
    dg3=handles.dg3;
    lenh=handles.lenh;
    
    set(handles.batch_button, 'Enable', 'off');
    
    if fnd>0

        for f=1:fnd
            
            axes(handles.progress2);
            cla;
            hold on;            
            bx=[0 f f 0 0]/fnd;
            by=[0 0 1 1 0];
            fill(bx,by,[.7 .7 .7]);            
            
            
            if f<=2 && fnd>3
                sdata{f,2}=' Reading images';
            elseif f>2 && size(sdata,1)>=3                
                sdata{3,2}=' Reading images';
            elseif fnd<=3
                sdata{f,2}=' Reading images';
            end                            
            
            set(handles.uitable1,'Data',sdata);
            pause(.001);

            okey=0;

            if strcmp(handles.ftype,'mat')
                bnd=load(efile{f},'boundaries');
                if ~isempty(bnd)
                    okey=1;
                    b=bnd.boundaries;
                    NUM_ELEMENTS=length(b);                
                    for k = 1:NUM_ELEMENTS               
                       [tmp,~]=regbound(b{k});
                       b{k}=tmp;                                                     
                    end                
                else
                    NUM_ELEMENTS=0;
                end

            elseif strcmp(handles.ftype,'tif')
                info = imfinfo(efile{f});
                NUM_ELEMENTS = numel(info);
                if NUM_ELEMENTS>0
                    okey=1;

                    b = cell(1,NUM_ELEMENTS);
                    for k = 1:NUM_ELEMENTS               
                        image = imread(efile{f},k);
                        if length(size(image))==3
                            image=image(:,:,1);
                        end
                        image = cast(image,'double');
                        image(image>0) = 1;
                        if max(image(:))==min(image(:))
                            image=ones(size(image));
                        end
                        [B,~]=bwboundaries(image,'noholes');
                        A=regionprops(logical(image),'Area');
                        Area=[A.Area];
                        [~,ind]=max(Area);
                        if ~isempty(ind)
                            [tmp,~]=regbound(B{ind});                        
                            b{k}=tmp;                                                                                              
                        end

                    end
                end            
            end

            if okey 
                
                boundaries=cell(1,NUM_ELEMENTS); sXbYb=boundaries; 
                sREM1=boundaries; sREM2=boundaries; sNconvex=boundaries; 
                sNconcave=boundaries; sNedge=boundaries; 
                sINDX=boundaries; sVb=boundaries; sCoor=boundaries; sDIST=boundaries; 
                sRad=boundaries; sTip=boundaries; sPathInd=boundaries; sPathLng=boundaries;
                
                for fr=1:NUM_ELEMENTS
                    
                    axes(handles.progress1);
                    cla;
                    hold on;                    
                    bx=[0 fr fr 0 0]/NUM_ELEMENTS;
                    by=[0 0 1 1 0];
                    fill(bx,by,[.7 .7 .7]);
                    pause(0.001);     
                     
                    if f<=2 && fnd>3
                        sdata{f,2}=[' Processing Frame ' num2str(fr)];
                    elseif f>2 && size(sdata,1)>=3                
                        sdata{3,2}=[' Processing Frame ' num2str(fr)];
                    elseif fnd<=3
                        sdata{f,2}=[' Processing Frame ' num2str(fr)];
                    end     
                                        
                    set(handles.uitable1,'Data',sdata);
                    pause(.001);

                    bound=round(b{fr});
                    by=bound(:,1); bx=bound(:,2); 

                    Status=0; atmp=0; tm=0;

                    while ~Status && atmp<30
                        atmp=atmp+1;       

                        Xb=bx+0.02*(randn(size(bx))-0.5);
                        Yb=by+0.02*(randn(size(bx))-0.5);

                        try
                            tic
                            [REM1,REM2,Status,Nconvex,Nconcave,Nedge]=FindGraphY(Xb',Yb',length(bx),handles);
                            tm=tm+toc;
                        catch
                            Status=0;
                        end
                    end
                    disp(' ');
                    disp(['Frame # ' num2str(fr) '; Attempts ' num2str(atmp) '; took ' num2str(round(tm)) ' sec.']);

                    if Status
                        REM=[REM1; REM2];
                        Smoo=handles.GN;
                        [INDX,Vb,Coor,DIST1,DIST2,Rad,~,~,Tip]=FindProfile(REM,Xb,Yb,Nconvex,Nconcave,Nedge,Smoo);
                        [PathInd,PathLng]=PathIndex(INDX,Vb);
                        DIST=(DIST1+DIST2)/2;

                        Lt=length(Tip);
                        tuck=zeros(1,Lt);
                        for i=1:Lt
                            [~,filo,~]=FindPathCond(INDX,Coor,DIST,Rad,Vb(Tip(i)),handles.CrR,handles.CutOff);
                            if filo(end,1)>0
                                tuck(i)=filo(end,1);
                            else
                                tuck(i)=filo(end-1,1);
                            end
                        end

                        nXb=zeros(size(Xb));
                        nYb=nXb; nLb=0;
                        for i=1:length(Xb)
                            tpath=PathInd(1:PathLng(i),i);
                            for j=1:Lt
                                tmp=find(tpath==tuck(j),1);
                                if ~isempty(tmp)
                                    break;
                                end
                            end
                            if isempty(tmp)
                                nLb=nLb+1;
                                nXb(nLb)=Xb(i);
                                nYb(nLb)=Yb(i);
                            end
                        end
                          
                        boundaries{fr}=[[nYb(1:nLb); nYb(1)],[nXb(1:nLb); nXb(1)]];
                        sXbYb{fr}=[Yb,Xb]; sREM1{fr}=REM1; sREM2{fr}=REM2; 
                        sNconvex{fr}=Nconvex; sNconcave{fr}=Nconcave; sNedge{fr}=Nedge; 
                        sINDX{fr}=INDX; sVb{fr}=Vb; sCoor{fr}=Coor; sDIST{fr}=DIST; sRad{fr}=Rad; sTip{fr}=Tip;
                        sPathInd{fr}=PathInd; sPathLng{fr}=PathLng;
                    else
                        disp(' ');
                        disp('Status = 0');
                    end

                end
                                
                if dg3>0
                    ins=efile{f}((end-3-dg3):(end-4));                
                    save([efile{f}(1:(end-lenh)) 'BGresult' ins '.mat'],...
                        'boundaries','sXbYb','sREM1','sREM2','sNconvex','sNconcave','sNedge',...
                        'sINDX','sVb','sCoor','sDIST','sRad','sTip','sPathInd','sPathLng');                
                else
                    save([efile{f}(1:(end-lenh)) 'BGresult.mat'],...
                        'boundaries','sXbYb','sREM1','sREM2','sNconvex','sNconcave','sNedge',...
                        'sINDX','sVb','sCoor','sDIST','sRad','sTip','sPathInd','sPathLng');                
                end
                
                
                if f<=2 && fnd>3
                    sdata{f,2}='   Processed';          
                elseif f>2 && size(sdata,1)>3
                    sdata{3,2}='   Processed';
                    sdata=sdata(2:end,:);                    
                elseif fnd>3 && size(sdata,1)==3
                    sdata{3,2}='   Processed';
                elseif fnd<=3
                    sdata{f,2}='   Processed';                    
                end    
                set(handles.uitable1,'Data',sdata);


            end
        end
        
        cla(handles.progress1);
        set(handles.progress1,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[],'Visible','off');
        cla(handles.progress2);
        set(handles.progress2,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[],'Visible','off');
        
        sdata=cell(1,2);
        sdata{1,1}='   All Files and Frames Are';       
        sdata{1,2}='   Processed';                    
        set(handles.uitable1,'Data',sdata);
        
    end
    
end
set(handle(jWindow),'Enabled',true);
%--------------------------------------------------------------------------   



% --- Checkboxes ----------------------------------------------------------
function checkbox1_Callback(hObject, eventdata, handles)

jFigPeer = get(handle(gcf),'JavaFrame');
try
    jWindow = jFigPeer.fFigureClient.getWindow;
catch
    jWindow = jFigPeer.fHG1Client.getWindow;
end
set(handle(jWindow),'Enabled',false);
try
    
    set(handles.movieSlider,'Enable','off');
    set(handles.lagSlider,'Enable','off');
    set(handles.CrR_Slider,'Enable','off');
    set(hObject,'Enable','off');
    
    handles = guidata(hObject);     
    handles = replot(handles);
       
    if handles.NUM_ELEMENTS>1
        set(handles.movieSlider,'Enable','on');
    end
    set(handles.lagSlider,'Enable','on');
    if get(hObject,'Value')
        set(handles.CrR_Slider,'Enable','on');
    end
    set(hObject,'Enable','on');
    guidata(hObject, handles);

end
set(handle(jWindow),'Enabled',true);    
%--------------------------------------------------------------------------


% --- Editing -------------------------------------------------------------
function edit_smoo_Callback(hObject, eventdata, handles)

    smoo=round(abs(str2double(get(hObject,'String'))));
    fr = round(get(handles.movieSlider,'Value'));
    bnd=handles.BOUNDARIES{fr};
    if smoo<size(bnd,1)
        handles.GN=smoo;
    end
    set(hObject,'String',num2str(handles.GN));
    if handles.DonePR(fr)
        try
            delete(handles.aLINE);
        end
        try
            delete(handles.bLINE);
        end
        %if handles.cLINE~=0 && ishandle(handles.cLINE) 
        try
            delete(handles.cLINE);            
        end
        axes(handles.profile_axis);        
        cla; axis fill;    
        set(handles.profile_axis,'Box','on','XTick',[],'YTick',[],'Visible','on'); 
        
        axes(handles.length_width_axis);        
        cla; axis fill;    
        set(handles.length_width_axis,'Box','on','XTick',[],'YTick',[],'Visible','on'); 
    end
    
    handles.DonePR = zeros(1,handles.NUM_ELEMENTS);
    set(handles.profile_button, 'Enable', 'on');
    set(handles.checkbox1, 'Enable', 'off');
    set(handles.lagSlider, 'Enable', 'off');
    set(handles.CrR_Slider, 'Enable', 'off');
    set(handles.CrR_text,'String','  ');
    set(handles.Tip_text,'String','  ');
    guidata(hObject, handles);
    %profile_button_Callback(hObject, eventdata, handles)
    
function edit_smoo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_smoo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% -------------------------------------------------------------------------



% --- Menu ----------------------------------------------------------------
function file_Callback(hObject, eventdata, handles)

function import_data_Callback(hObject, eventdata, handles)

    import(hObject,2);
%{
function import_mat_Callback(hObject, eventdata, handles)

    import(hObject,1);

function import_tif_Callback(hObject, eventdata, handles)

    import(hObject,2);
%}   
function save_res_Callback(hObject, eventdata, handles)

    fr = round(get(handles.movieSlider,'Value'));
    if sum(handles.DoneGR)==handles.NUM_ELEMENTS && sum(handles.DonePR)==handles.NUM_ELEMENTS
        
        boundaries=handles.gbound;
        sXbYb=handles.BOUNDARIES;
        sREM1=handles.REM1;
        sREM2=handles.REM2;
        sNconvex=handles.Nconvex;
        sNconcave=handles.Nconcave;
        sNedge=handles.Nedge;
        
        sINDX=handles.INDX;
        sVb=handles.Vb;
        sCoor=handles.Coor;
        sDIST=handles.DIST;
        sRad=handles.Rad;
        sTip=handles.Tip;
            
        sPathInd=handles.PathInd;
        sPathLng=handles.PathLng;    
        
        uisave({'boundaries','sXbYb','sREM1','sREM2','sNconvex','sNconcave','sNedge',...
                'sINDX','sVb','sCoor','sDIST','sRad','sTip','sPathInd','sPathLng'},'WholeMovieData');            
        
        
    elseif handles.DoneGR(fr) && handles.DonePR(fr)
        boundaries=handles.gbound(fr);
        sXbYb=handles.BOUNDARIES(fr);
        sREM1=handles.REM1(fr);
        sREM2=handles.REM2(fr);
        sNconvex=handles.Nconvex(fr);
        sNconcave=handles.Nconcave(fr);
        sNedge=handles.Nedge(fr);
        
        sINDX=handles.INDX{fr};
        sVb=handles.Vb{fr};
        sCoor=handles.Coor{fr};
        sDIST=handles.DIST{fr};
        sRad=handles.Rad{fr};
        sTip=handles.Tip{fr};
            
        sPathInd=handles.PathInd{fr};
        sPathLng=handles.PathLng{fr};    
        
        uisave({'boundaries','sXbYb','sREM1','sREM2','sNconvex','sNconcave','sNedge',...
                'sINDX','sVb','sCoor','sDIST','sRad','sTip','sPathInd','sPathLng'},'CurrentFrameData');
    end
% -------------------------------------------------------------------------
    


function lev_1_1_Callback(hObject, eventdata, handles)

function lev_1_2_Callback(hObject, eventdata, handles)

function lev_1_3_Callback(hObject, eventdata, handles)

function lev_1_4_Callback(hObject, eventdata, handles)

function lev_1_5_Callback(hObject, eventdata, handles)

function lev_2_1_Callback(hObject, eventdata, handles)

function lev_2_2_Callback(hObject, eventdata, handles)

function lev_2_3_Callback(hObject, eventdata, handles)

function lev_2_4_Callback(hObject, eventdata, handles)

function lev_2_5_Callback(hObject, eventdata, handles)

function lev_3_1_Callback(hObject, eventdata, handles)

function lev_3_2_Callback(hObject, eventdata, handles)

function lev_3_3_Callback(hObject, eventdata, handles)

function lev_3_4_Callback(hObject, eventdata, handles)

function lev_3_5_Callback(hObject, eventdata, handles)


function lev_1_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lev_1_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lev_1_3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lev_1_4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lev_1_5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lev_2_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lev_2_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lev_2_3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lev_2_4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lev_2_5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lev_3_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lev_3_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lev_3_3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lev_3_4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lev_3_5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
