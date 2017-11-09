function varargout = ProActive(varargin)
% PROACTIVE MATLAB code for ProActive.fig
%      PROACTIVE, by itself, creates a new PROACTIVE or raises the existing
%      singleton*.
%
%      H = PROACTIVE returns the handle to a new PROACTIVE or the handle to
%      the existing singleton*.
%
%      PROACTIVE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROACTIVE.M with the given input arguments.
%
%      PROACTIVE('Property','Value',...) creates a new PROACTIVE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ProActive_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ProActive_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ProActive

% Last Modified by GUIDE v2.5 12-Feb-2013 10:21:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ProActive_OpeningFcn, ...
                   'gui_OutputFcn',  @ProActive_OutputFcn, ...
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

function ProActive_OpeningFcn(hObject, eventdata, handles, varargin)

    handles.output = hObject;
    def=importdata('Defaults.dat');

    handles.NUM_ELEMENTS = 0;
    handles.IMAGE_NUMBER = 1;
    handles.X = 0; 
    handles.Y = 0; 
    
    handles.Lag = 1;
    handles.NeedSpeed = 0;
    handles.DefFilt = def(1);
    handles.DefRange = def(2);
    handles.DefRange2 = def(3);
    handles.DefGN = def(4);
    
    handles.MeanVelPer=0;
    handles.MaxVelPer=0;
    handles.zMaxVelPer=0;
    handles.TVelPer=0;
    handles.zTVelPer=0;

    % Listeners
    handles.LisMov = addlistener(handles.movieSlider, 'Value', 'PostSet', @(src, event)switchFrame(hObject, src, event));
    handles.LisLag = addlistener(handles.lagSlider, 'Value', 'PostSet', @(src, event)switchLag(hObject, src, event));
    handles.LisFil = addlistener(handles.filterSlider, 'Value', 'PostSet', @(src, event)switchFilt(hObject, src, event));

    set(handles.frame_text,'ForegroundColor',[.5 .5 .5]);
    set(handles.lag_text,'ForegroundColor',[.5 .5 .5]);
    set(handles.filter_text,'ForegroundColor',[.5 .5 .5]);
    set(handles.mean_text,'ForegroundColor',[.5 .5 .5]);
    set(handles.max_text,'ForegroundColor',[.5 .5 .5]);
    set(handles.peak_text,'ForegroundColor',[.5 .5 .5]);
    set(handles.upto1,'ForegroundColor',[.5 .5 .5]);
    set(handles.upto2,'ForegroundColor',[.5 .5 .5]);
    set(handles.activity_type,'ForegroundColor',[.5 .5 .5]);
    set(handles.norm_type,'ForegroundColor',[.5 .5 .5]);
    set(handles.vel_panel,'ForegroundColor',[.5 .5 .5]);
    set(handles.results_panel,'ForegroundColor',[.5 .5 .5]);
    set(handles.smoo_text,'ForegroundColor',[.5 .5 .5]);

    % Update handles structure
    set(handles.boundary_axis,'Box','on','XTick',[],'YTick',[]);
    set(handles.activity_axis,'Box','on','XTick',[],'YTick',[]);
    set(handles.results_axis,'Box','on','XTick',[],'YTick',[]);

    set(handles.movieSlider, 'Enable', 'off');
    set(handles.lagSlider, 'Enable', 'off');
    set(handles.filterSlider, 'Enable', 'off');
    set(handles.run_lag_range, 'Enable', 'off');

    set(handles.pro_act, 'Enable', 'off');
    set(handles.retro_act, 'Enable', 'off');
    set(handles.total_act, 'Enable', 'off');
    set(handles.perim_norm, 'Enable', 'off');
    set(handles.area_norm, 'Enable', 'off');
    set(handles.none_norm, 'Enable', 'off');

    set(handles.checkbox1, 'Enable', 'off');
    set(handles.checkbox2, 'Enable', 'off');

    set(handles.save_as, 'Enable', 'off');

    guidata(hObject, handles);

function varargout = ProActive_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function main_window_CloseRequestFcn(hObject, eventdata, handles)

    try 
        def = [handles.Filt, handles.Range, handles.Range2, handles.GN];
        id = fopen([cd '/Defaults.dat'],'w');
        fprintf(id,'%d \n',def);  
        close(id);
    end
    delete(hObject);    
%--------------------------------------------------------------------------



% --- Function for Listeners ---------------------------------------------- 
function switchFrame(hObject, src, event)

    handles = guidata(hObject);    
    sliderNumber = round(event.NewValue);    
    handles.IMAGE_NUMBER = sliderNumber;
    
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
        
    axes(handles.activity_axis);
    xlim=get(handles.activity_axis,'XLim');
    ylim=get(handles.activity_axis,'YLim');
    cla;
    hold on;    
    if (sliderNumber+handles.Lag)<=handles.NUM_ELEMENTS
        M=handles.mask{sliderNumber+handles.Lag}-handles.mask{sliderNumber};
    else
        M=handles.empty;
    end
    imagesc(M, 'Parent', handles.activity_axis);
    %{
            bound=handles.BOUNDARIES{sliderNumber};
            by=bound(:,1); bx=bound(:,2);
            plot(bx-handles.X(1)+1,by-handles.Y(1)+1,'m');
            bound=handles.BOUNDARIES{sliderNumber+handles.Lag};
            by=bound(:,1); bx=bound(:,2);
            plot(bx-handles.X(1)+1,by-handles.Y(1)+1,'c');
    %}
    if handles.NeedSpeed && (sliderNumber+handles.Lag)<=handles.NUM_ELEMENTS
        bound=handles.BOUNDARIES{sliderNumber+handles.Lag};
        by=bound(:,1); bx=bound(:,2);
        ind=(handles.Speed{sliderNumber,handles.Lag}>handles.Filt);        
        plot(bx(ind)-handles.X(1)+1,by(ind)-handles.Y(1)+1,'r.');
    end
    axis image;
    axis ij;
    set(handles.activity_axis,'XLim',xlim);
    set(handles.activity_axis,'YLim',ylim);
    axis off;
   
    set(handles.frame_text,'String',['Frame # ' num2str(sliderNumber)]);
    
    guidata(hObject, handles);
   
function switchLag(hObject, src, event) 
    handles = guidata(hObject);
    lag = round(event.NewValue);
    handles.Lag = lag;
    set(handles.lag_text, 'String', ['Lag = ' num2str(lag)]);
    guidata(hObject, handles);
   
function switchFilt(hObject, src, event) 
    handles = guidata(hObject);
    filt = round(event.NewValue);
    handles.Filt = filt;
    set(handles.filter_text, 'String', ['Filter = ' num2str(filt)]);
    
    if handles.Done2(handles.Lag)
    
        axes(handles.activity_axis);
        xlim=get(handles.activity_axis,'XLim');
        ylim=get(handles.activity_axis,'YLim');
        cla;
        hold on;    
        if (handles.IMAGE_NUMBER+handles.Lag)<=handles.NUM_ELEMENTS
            M=handles.mask{handles.IMAGE_NUMBER+handles.Lag}-handles.mask{handles.IMAGE_NUMBER};
        else
            M=handles.empty;
        end
        imagesc(M,'Parent', handles.activity_axis);

        if handles.NeedSpeed && (handles.IMAGE_NUMBER+handles.Lag)<=handles.NUM_ELEMENTS
            bound=handles.BOUNDARIES{handles.IMAGE_NUMBER+handles.Lag};
            by=bound(:,1); bx=bound(:,2);
            ind=(handles.Speed{handles.IMAGE_NUMBER,handles.Lag}>filt);
            plot(bx(ind)-handles.X(1)+1,by(ind)-handles.Y(1)+1,'r.');          
        end

        axis image;
        axis ij;
        set(handles.activity_axis,'XLim',xlim);
        set(handles.activity_axis,'YLim',ylim);
        axis off;
    end
    
    guidata(hObject, handles);
%--------------------------------------------------------------------------       



% --- Additional functions ------------------------------------------------
function import(hObject,imp)
    handles = guidata(hObject);
    
    handles.Lag = 1;
    handles.IMAGE_NUMBER = 1;
    handles.NeedSpeed = 0;
    delete(handles.LisMov);
    delete(handles.LisLag);
    delete(handles.LisFil);
    
    axes(handles.activity_axis);
    %set(handles.activity_axis,'Position',[550 295 500 430]);
    colormap(gray); 
    cla; axis fill;    
    set(handles.activity_axis,'Box','on','XTick',[],'YTick',[],'Visible','on');
    
    axes(handles.boundary_axis);
    %set(handles.boundary_axis,'Position',[20 295 500 430]);
    cla; axis fill;    
    set(handles.boundary_axis,'Box','on','XTick',[],'YTick',[],'Visible','on');
    
    axes(handles.results_axis);
    %set(handles.results_axis,'Position',[50 35 470 180]);
    cla; axis fill;    
    set(handles.results_axis,'Box','on','XTick',[],'YTick',[],'Visible','on');
    
    set(handles.lag_text, 'String', 'Lag','ForegroundColor',[.5 .5 .5]);
    set(handles.frame_text,'String','Frame','ForegroundColor',[.5 .5 .5]);
    set(handles.filter_text,'String','Filter','ForegroundColor',[.5 .5 .5]);
    set(handles.mean_text, 'String', 'Mean','ForegroundColor',[.5 .5 .5]);
    set(handles.max_text, 'String', 'Max','ForegroundColor',[.5 .5 .5]);
    set(handles.peak_text, 'String', 'at T ','ForegroundColor',[.5 .5 .5]);    
    set(handles.edit_range,'String',' ','Enable','off');
    set(handles.edit_range2,'String',' ','Enable','off');
    set(handles.edit_smoo,'String',' ','Enable','off');
    set(handles.upto1,'ForegroundColor',[.5 .5 .5]);
    set(handles.upto2,'ForegroundColor',[.5 .5 .5]);
    set(handles.activity_type,'ForegroundColor',[.5 .5 .5]);
    set(handles.norm_type,'ForegroundColor',[.5 .5 .5]);
    set(handles.vel_panel,'ForegroundColor',[.5 .5 .5]);
    set(handles.results_panel,'ForegroundColor',[.5 .5 .5]);
    set(handles.smoo_text,'ForegroundColor',[.5 .5 .5]);
    set(handles.result_title,'String','  ');      
    
    set(handles.movieSlider, 'Enable', 'off');
    set(handles.lagSlider, 'Enable', 'off');
    set(handles.filterSlider, 'Enable', 'off');
    set(handles.run_lag_range, 'Enable', 'off');
    
    set(handles.pro_act, 'Enable', 'off');
    set(handles.retro_act, 'Enable', 'off');
    set(handles.total_act, 'Enable', 'off');
    set(handles.perim_norm, 'Enable', 'off');
    set(handles.area_norm, 'Enable', 'off');
    set(handles.none_norm, 'Enable', 'off');
    
    set(handles.checkbox1, 'Value', 0);
    set(handles.checkbox1, 'Enable', 'off');
    set(handles.checkbox2, 'Value', 0);
    set(handles.checkbox2, 'Enable', 'off');
    
    set(handles.save_as, 'Enable', 'off');

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

                    mask = cell(1,handles.NUM_ELEMENTS);
                    TArea = zeros(1,handles.NUM_ELEMENTS);
                    TPerm = zeros(1,handles.NUM_ELEMENTS);
                    xmi= zeros(1,handles.NUM_ELEMENTS);
                    xma=xmi; ymi=xmi; yma=xmi;

                    for k = 1:handles.NUM_ELEMENTS               
                       tmp=round(b{k}); 
                       b{k}=tmp;               
                       xmi(k)=min(tmp(:,2)); xma(k)=max(tmp(:,2));
                       ymi(k)=min(tmp(:,1)); yma(k)=max(tmp(:,1));                
                    end

                    xlim=[min(xmi)-10 max(xma)+10];
                    ylim=[min(ymi)-10 max(yma)+10];
                    handles.empty = zeros(ylim(2)-ylim(1)+1,xlim(2)-xlim(1)+1);

                    set(handles.frame_text,'ForegroundColor','k');
                    for k = 1:handles.NUM_ELEMENTS
                       set(handles.frame_text,'String',['Frame # ' num2str(k)]); 
                       pause(.001); 
                       tmp=b{k};
                       p2m=poly2mask(tmp(:,2)-xlim(1)+1, tmp(:,1)-ylim(1)+1, ylim(2)-ylim(1)+1, xlim(2)-xlim(1)+1);
                       p2m((ylim(2)-ylim(1)+1)*(tmp(:,2)-xlim(1))+tmp(:,1)-ylim(1)+1)=1;
                       mask{k}=p2m;
                       TArea(k)=polyarea(tmp(:,2),tmp(:,1));               
                       tleng=sqrt((tmp(2:end,2)-tmp(1:(end-1),2)).^2+(tmp(2:end,1)-tmp(1:(end-1),1)).^2);
                       TPerm(k)=sum(tleng)+sqrt((tmp(1,2)-tmp(end,2))^2+(tmp(1,1)-tmp(end,1))^2);                              
                    end
                end

            elseif imp==2
                info = imfinfo([pathname filename]);
                handles.NUM_ELEMENTS = numel(info);
                if handles.NUM_ELEMENTS>1
                    okey=1;

                    b = cell(1,handles.NUM_ELEMENTS);
                    mask = cell(1,handles.NUM_ELEMENTS);
                    TArea = zeros(1,handles.NUM_ELEMENTS);
                    TPerm = zeros(1,handles.NUM_ELEMENTS);
                    xmi= zeros(1,handles.NUM_ELEMENTS);
                    xma=xmi; ymi=xmi; yma=xmi;

                    set(handles.frame_text,'ForegroundColor','k');
                    for k = 1:handles.NUM_ELEMENTS 
                        set(handles.frame_text,'String',['Frame # ' num2str(k)]); 
                        pause(.001);
                        image = imread([pathname filename], k);
                        if length(size(image))==3
                            image=image(:,:,1);
                        end
                        image = cast(image,'double');
                        image(image>0) = 1;
                        if max(image(:))==min(image(:))
                            image=ones(size(image));
                        end

                        [B,L]=bwboundaries(image,'noholes');
                        A=regionprops(logical(image),'Area');
                        Area=[A.Area];
                        [v,ind]=max(Area);
                        if ~isempty(ind)
                            tmp=B{ind}; L(L~=ind)=0; L(L==ind)=1;
                            smx=sum(L,1); smy=sum(L,2);
                            xmi(k)=find(smx>0,1,'first');
                            xma(k)=find(smx>0,1,'last');
                            ymi(k)=find(smy>0,1,'first');
                            yma(k)=find(smy>0,1,'last');
                            b{k}=tmp;
                            TArea(k)=v;
                            tleng=sqrt((tmp(2:end,2)-tmp(1:(end-1),2)).^2+(tmp(2:end,1)-tmp(1:(end-1),1)).^2);
                            TPerm(k)=sum(tleng)+sqrt((tmp(1,2)-tmp(end,2))^2+(tmp(1,1)-tmp(end,1))^2);                                                      
                        else
                            xmi(k)=Inf;
                            xma(k)=-Inf;
                            ymi(k)=Inf;
                            yma(k)=-Inf;
                        end
                        mask{k}=imfill(L,'holes');                    
                    end

                    xlim=[min(xmi)-10 max(xma)+10];
                    ylim=[min(ymi)-10 max(yma)+10];
                    handles.empty = zeros(ylim(2)-ylim(1)+1,xlim(2)-xlim(1)+1);


                    for k = 1:handles.NUM_ELEMENTS                                            
                        imagec=handles.empty;                    
                        imagec(11:(ylim(2)-ylim(1)-9),11:(xlim(2)-xlim(1)-9))...
                            = mask{k}((ylim(1)+10):(ylim(2)-10),(xlim(1)+10):(xlim(2)-10));
                        mask{k} = imagec;                   
                    end                
                else
                    errordlg('The data consists of 1 time frame. The program is meant to analysis movies, not single frames.');
                end            
            end

            if okey

                if handles.DefRange < handles.NUM_ELEMENTS
                    handles.Range = handles.DefRange;
                else
                    handles.Range = handles.NUM_ELEMENTS-1;
                end
                if handles.DefGN < handles.NUM_ELEMENTS
                    handles.GN = handles.DefGN;
                else
                    handles.GN = handles.NUM_ELEMENTS-1;
                end            
                handles.Filt = handles.DefFilt;
                handles.Range2 =  handles.DefRange2;

                set(handles.movieSlider, 'Min', 1, 'Max', handles.NUM_ELEMENTS,'Value',1);
                if handles.NUM_ELEMENTS>2
                    set(handles.lagSlider,'Min', 1, 'Max', handles.NUM_ELEMENTS-1,'Value',1);
                end
                set(handles.filterSlider,'Min', 0, 'Max', 50,'Value',handles.Filt);

                handles.PrActF = cell(1,handles.NUM_ELEMENTS-1);
                handles.PrActB = cell(1,handles.NUM_ELEMENTS-1);
                handles.PrActT = cell(1,handles.NUM_ELEMENTS-1);
                handles.PrActR = cell(1,handles.NUM_ELEMENTS-1);
                handles.Done1 = zeros(1,handles.NUM_ELEMENTS-1);
                handles.Done2 = zeros(1,handles.NUM_ELEMENTS-1);

                PrActF = zeros(1,handles.NUM_ELEMENTS-handles.Lag);
                PrActB = zeros(1,handles.NUM_ELEMENTS-handles.Lag);
                PrActT = zeros(1,handles.NUM_ELEMENTS-handles.Lag);
                PrActR = zeros(6,handles.NUM_ELEMENTS-handles.Lag);

                range=handles.Range;
                handles.MeanActF=zeros(1,range); handles.MeanActB=zeros(1,range); handles.MeanActT=zeros(1,range); handles.MeanActR=zeros(6,range);
                handles.MaxActF=zeros(1,range);  handles.MaxActB=zeros(1,range);  handles.MaxActT=zeros(1,range);  handles.MaxActR=zeros(6,range);
                handles.TActF=zeros(1,range);    handles.TActB=zeros(1,range);    handles.TActT=zeros(1,range);    handles.TActR=zeros(6,range);
                handles.zMaxActF=zeros(1,range); handles.zMaxActB=zeros(1,range); handles.zMaxActT=zeros(1,range); handles.zMaxActR=zeros(6,range);
                handles.zTActF=zeros(1,range);   handles.zTActB=zeros(1,range);   handles.zTActT=zeros(1,range);   handles.zTActR=zeros(6,range);

                for k = 1:(handles.NUM_ELEMENTS-handles.Lag)                              
                    M=mask{k+handles.Lag}-mask{k};                
                    PrActF(k)=sum(sum(M==+1));
                    PrActB(k)=sum(sum(M==-1));
                    PrActT(k)=PrActF(k)+PrActB(k); 
                    PrActR(1,k)=100*2*PrActF(k)/(TArea(k)+TArea(k+handles.Lag));
                    PrActR(2,k)=100*2*PrActB(k)/(TArea(k)+TArea(k+handles.Lag));
                    PrActR(3,k)=100*2*PrActT(k)/(TArea(k)+TArea(k+handles.Lag));
                    PrActR(4,k)=    2*PrActF(k)/(TPerm(k)+TPerm(k+handles.Lag));
                    PrActR(5,k)=    2*PrActB(k)/(TPerm(k)+TPerm(k+handles.Lag));
                    PrActR(6,k)=    2*PrActT(k)/(TPerm(k)+TPerm(k+handles.Lag));                
                end           

                handles.mask = mask;
                handles.BOUNDARIES = b;
                handles.X = xlim;
                handles.Y = ylim;

                handles.TArea=TArea;
                handles.TPerm=TPerm;
                handles.PrActF{handles.Lag}=PrActF;
                handles.PrActB{handles.Lag}=PrActB;
                handles.PrActT{handles.Lag}=PrActT; 
                handles.PrActR{handles.Lag}=PrActR; 

                zPrActF=GFilterA(PrActF,handles.GN); zPrActB=GFilterA(PrActB,handles.GN); 
                zPrActT=GFilterA(PrActT,handles.GN); zPrActR=PrActR;
                zPrActR(1,:)=GFilterA(PrActR(1,:),handles.GN); zPrActR(2,:)=GFilterA(PrActR(2,:),handles.GN);
                zPrActR(3,:)=GFilterA(PrActR(3,:),handles.GN); zPrActR(4,:)=GFilterA(PrActR(4,:),handles.GN);
                zPrActR(5,:)=GFilterA(PrActR(5,:),handles.GN); zPrActR(6,:)=GFilterA(PrActR(6,:),handles.GN);

                handles.MeanActF(handles.Lag)=mean(PrActF);
                handles.MeanActB(handles.Lag)=mean(PrActB);
                handles.MeanActT(handles.Lag)=mean(PrActT);
                handles.MeanActR(:,handles.Lag)=mean(PrActR,2);

                [handles.MaxActF(handles.Lag),handles.TActF(handles.Lag)]=max(PrActF);
                [handles.MaxActB(handles.Lag),handles.TActB(handles.Lag)]=max(PrActB);
                [handles.MaxActT(handles.Lag),handles.TActT(handles.Lag)]=max(PrActT);
                [handles.MaxActR(:,handles.Lag),handles.TActR(:,handles.Lag)]=max(PrActR,[],2);  

                [handles.zMaxActF(handles.Lag),handles.zTActF(handles.Lag)]=max(zPrActF);
                [handles.zMaxActB(handles.Lag),handles.zTActB(handles.Lag)]=max(zPrActB);
                [handles.zMaxActT(handles.Lag),handles.zTActT(handles.Lag)]=max(zPrActT);
                [handles.zMaxActR(:,handles.Lag),handles.zTActR(:,handles.Lag)]=max(zPrActR,[],2);  

                handles.Done1(handles.Lag)=1;

                if(get(handles.pro_act,'Value') == 1)
                    if(get(handles.perim_norm, 'Value') == 1)            
                        Y = PrActR(4,:);       
                        tit='Protrusive area normalized by cell perimeter';
                    elseif (get(handles.area_norm, 'Value') == 1)    
                        Y = PrActR(1,:);
                        tit='Protrusive area normalized by cell area';
                    elseif (get(handles.none_norm, 'Value') == 1)
                        Y = PrActF;
                        tit='Protrusive area (not normalized)';
                    end
                elseif(get(handles.retro_act,'Value') == 1)
                    if(get(handles.perim_norm, 'Value') == 1)
                        Y = PrActR(5,:);
                        tit='Retractive area normalized by cell perimeter';
                    elseif (get(handles.area_norm, 'Value') == 1)    
                        Y = PrActR(2,:);
                        tit='Retractive area normalized by cell area';
                    elseif (get(handles.none_norm, 'Value') == 1)
                        Y = PrActB;
                        tit='Retractive area (not normalized)';
                    end  
                elseif(get(handles.total_act,'Value') == 1)
                    if(get(handles.perim_norm, 'Value') == 1)
                        Y = PrActR(6,:);      
                        tit='Total area normalized by cell perimeter';
                    elseif (get(handles.area_norm, 'Value') == 1)    
                        Y = PrActR(3,:);
                        tit='Total area normalized by cell area';
                    elseif (get(handles.none_norm, 'Value') == 1)
                        Y = PrActT;
                        tit='Total area (not normalized)';
                    end
                end

                axes(handles.activity_axis);
                hold on;
                imagesc(handles.mask{2}-handles.mask{1});
                %{
                    bound=b{1};
                    by=bound(:,1); bx=bound(:,2);
                    plot(bx-xlim(1)+1,by-ylim(1)+1,'m');
                    bound=b{2};
                    by=bound(:,1); bx=bound(:,2);
                    plot(bx-xlim(1)+1,by-ylim(1)+1,'c');
                %}
                axis image; 
                axis ij;
                axis off;

                axes(handles.boundary_axis);
                hold on;
                bound=b{1};
                by=bound(:,1); bx=bound(:,2);
                fill(bx,by,[.9 .9 .9],'EdgeColor','k');  
                axis image; 
                axis ij;
                set(handles.boundary_axis,'XLim',xlim,'YLim',ylim);
                axis off;

                axes(handles.results_axis);         
                cla reset;
                hold on;            
                [mv mi]=max(Y);
                plot([1 handles.NUM_ELEMENTS],[mv mv],':','Color',[.6 .6 .6]);
                plot([mi mi],[0 mv],':','Color',[.6 .6 .6]);
                plot([1 handles.NUM_ELEMENTS],[mean(Y) mean(Y)],'r','LineWidth',2);
                plot(Y,'g','LineWidth',2);            
                if min(Y)<mv
                    set(handles.results_axis,'Box','on','XLim',[1 handles.NUM_ELEMENTS],'YLim',[1.1*min(Y)-0.1*mv 1.1*mv-0.1*min(Y)]);
                elseif mv>0
                    set(handles.results_axis,'Box','on','XLim',[1 handles.NUM_ELEMENTS],'YLim',[0.9*mv 1.1*mv]);
                else
                    set(handles.results_axis,'Box','on','XLim',[1 handles.NUM_ELEMENTS],'YLim',[-1 1]);
                end

                set(handles.mean_text,'String', ['Mean = ' num2str(mean(Y),'%.2f')],'ForegroundColor','k');
                set(handles.max_text, 'String', ['Max = ' num2str(mv,'%.2f')],'ForegroundColor','k');
                set(handles.peak_text,'String', ['at T = ' num2str(mi)],'ForegroundColor','k');   

                set(handles.result_title,'String',tit);
                set(handles.lag_text, 'String', 'Lag = 1','ForegroundColor','k');
                set(handles.frame_text,'String','Frame # 1');
                set(handles.edit_range,'String',num2str(handles.Range),'Enable','on');
                set(handles.edit_range2,'String',' ');

                set(handles.upto1,'ForegroundColor','k');
                set(handles.activity_type,'ForegroundColor','k');
                set(handles.norm_type,'ForegroundColor','k');
                set(handles.vel_panel,'ForegroundColor','k');
                set(handles.results_panel,'ForegroundColor','k');

                set(handles.movieSlider, 'Enable', 'on');
                if handles.NUM_ELEMENTS>2
                    set(handles.lagSlider, 'Enable', 'on');
                end
                set(handles.run_lag_range, 'Enable', 'on');

                set(handles.pro_act, 'Enable', 'on');
                set(handles.retro_act, 'Enable', 'on');
                set(handles.total_act, 'Enable', 'on');
                set(handles.perim_norm, 'Enable', 'on');
                set(handles.area_norm, 'Enable', 'on');
                set(handles.none_norm, 'Enable', 'on');    

                set(handles.checkbox1, 'Enable', 'on');
                set(handles.checkbox2, 'Enable', 'on');

                set(handles.save_as, 'Enable', 'on');

                guidata(hObject, handles);
            end
        catch
            errordlg('Import failed. Please check the file format.');    
        end
        set(handle(jWindow),'Enabled',true);
    end
    handles.LisMov = addlistener(handles.movieSlider, 'Value', 'PostSet', @(src, event)switchFrame(hObject, src, event));
    handles.LisLag = addlistener(handles.lagSlider, 'Value', 'PostSet', @(src, event)switchLag(hObject, src, event));
    handles.LisFil = addlistener(handles.filterSlider, 'Value', 'PostSet', @(src, event)switchFilt(hObject, src, event));

    guidata(hObject,handles);

function handles=LAG_STAT(handles)
           
    b=handles.BOUNDARIES;
    xlim=handles.X;
    ylim=handles.Y;
    TArea=handles.TArea;
    TPerm=handles.TPerm;
    
    PrActF = zeros(1,handles.NUM_ELEMENTS-handles.Lag);
    PrActB = zeros(1,handles.NUM_ELEMENTS-handles.Lag);
    PrActT = zeros(1,handles.NUM_ELEMENTS-handles.Lag);
    PrActR = zeros(6,handles.NUM_ELEMENTS-handles.Lag);
    
    set(handles.movieSlider, 'Enable', 'off');
    set(handles.lagSlider, 'Enable', 'off');
    for k = 1:(handles.NUM_ELEMENTS-handles.Lag)
        tmp=b{k};
        Xb1=tmp(:,2)-xlim(1)+1; Yb1=tmp(:,1)-ylim(1)+1;
        tmp=b{k+handles.Lag};
        Xb2=tmp(:,2)-xlim(1)+1; Yb2=tmp(:,1)-ylim(1)+1;
        M=handles.mask{k+handles.Lag}-handles.mask{k};
        if handles.NeedSpeed                       
            tmp=handles.mask{k};
            dist=zeros(1,length(Xb2));
            in=zeros(length(Xb2),1);
            for d=1:length(Xb2)
                dist(d)=min(abs(Xb2(d)-Xb1)+abs(Yb2(d)-Yb1));   
                %in(d)=1-2*tmp(round(Yb2(d)),round(Xb2(d)));
                in(d)=1-2*tmp(Yb2(d),Xb2(d));
            end
            handles.Speed{k,handles.Lag}=in'.*dist;            
        end
        PrActF(k)=sum(sum(M==+1));
        PrActB(k)=sum(sum(M==-1));
        PrActT(k)=PrActF(k)+PrActB(k); 
        PrActR(1,k)=100*2*PrActF(k)/(TArea(k)+TArea(k+handles.Lag));
        PrActR(2,k)=100*2*PrActB(k)/(TArea(k)+TArea(k+handles.Lag));
        PrActR(3,k)=100*2*PrActT(k)/(TArea(k)+TArea(k+handles.Lag));
        PrActR(4,k)=    2*PrActF(k)/(TPerm(k)+TPerm(k+handles.Lag));
        PrActR(5,k)=    2*PrActB(k)/(TPerm(k)+TPerm(k+handles.Lag));
        PrActR(6,k)=    2*PrActT(k)/(TPerm(k)+TPerm(k+handles.Lag));        
    end

    handles.PrActF{handles.Lag}=PrActF;
    handles.PrActB{handles.Lag}=PrActB;
    handles.PrActT{handles.Lag}=PrActT; 
    handles.PrActR{handles.Lag}=PrActR;
    
    zPrActF=GFilterA(PrActF,handles.GN); zPrActB=GFilterA(PrActB,handles.GN); 
    zPrActT=GFilterA(PrActT,handles.GN); zPrActR=PrActR;
    zPrActR(1,:)=GFilterA(PrActR(1,:),handles.GN); zPrActR(2,:)=GFilterA(PrActR(2,:),handles.GN);
    zPrActR(3,:)=GFilterA(PrActR(3,:),handles.GN); zPrActR(4,:)=GFilterA(PrActR(4,:),handles.GN);
    zPrActR(5,:)=GFilterA(PrActR(5,:),handles.GN); zPrActR(6,:)=GFilterA(PrActR(6,:),handles.GN);
    
    handles.MeanActF(handles.Lag)=mean(PrActF);
    handles.MeanActB(handles.Lag)=mean(PrActB);
    handles.MeanActT(handles.Lag)=mean(PrActT);
    handles.MeanActR(:,handles.Lag)=mean(PrActR,2);

    [handles.MaxActF(handles.Lag),handles.TActF(handles.Lag)]=max(PrActF);
    [handles.MaxActB(handles.Lag),handles.TActB(handles.Lag)]=max(PrActB);
    [handles.MaxActT(handles.Lag),handles.TActT(handles.Lag)]=max(PrActT);
    [handles.MaxActR(:,handles.Lag),handles.TActR(:,handles.Lag)]=max(PrActR,[],2);
    
    [handles.zMaxActF(handles.Lag),handles.zTActF(handles.Lag)]=max(zPrActF);
    [handles.zMaxActB(handles.Lag),handles.zTActB(handles.Lag)]=max(zPrActB);
    [handles.zMaxActT(handles.Lag),handles.zTActT(handles.Lag)]=max(zPrActT);
    [handles.zMaxActR(:,handles.Lag),handles.zTActR(:,handles.Lag)]=max(zPrActR,[],2);
    
    handles.Done1(handles.Lag)=1;
    if handles.NeedSpeed        
        handles.Done2(handles.Lag)=1;
    end
    set(handles.movieSlider, 'Enable', 'on');
    set(handles.lagSlider, 'Enable', 'on');
%--------------------------------------------------------------------------    
  


% --- Sliders -------------------------------------------------------------
function movieSlider_Callback(hObject, eventdata, handles)

function lagSlider_Callback(hObject, eventdata, handles)

    lag = round(get(handles.lagSlider,'Value'));
    
    jFigPeer = get(handle(gcf),'JavaFrame');
    try
        jWindow = jFigPeer.fFigureClient.getWindow;
    catch
        jWindow = jFigPeer.fHG1Client.getWindow;
    end
    set(handle(jWindow),'Enabled',false);
    try
        if handles.NeedSpeed && handles.Done2(lag)==0
            handles=LAG_STAT(handles);        
        elseif ~handles.NeedSpeed && handles.Done1(lag)==0
            handles=LAG_STAT(handles);               
        end

        axes(handles.activity_axis);
        xlim=get(handles.activity_axis,'XLim');
        ylim=get(handles.activity_axis,'YLim');
        cla;
        hold on;    
        if (handles.IMAGE_NUMBER+lag)<=handles.NUM_ELEMENTS
            M=handles.mask{handles.IMAGE_NUMBER+lag}-handles.mask{handles.IMAGE_NUMBER};
        else
            M=handles.empty;
        end
        imagesc(M,'Parent', handles.activity_axis);

        if handles.NeedSpeed && (handles.IMAGE_NUMBER+lag)<=handles.NUM_ELEMENTS
            bound=handles.BOUNDARIES{handles.IMAGE_NUMBER+lag};
            by=bound(:,1); bx=bound(:,2);
            ind=(handles.Speed{handles.IMAGE_NUMBER,lag}>handles.Filt);
            plot(bx(ind)-handles.X(1)+1,by(ind)-handles.Y(1)+1,'r.');          
        end

        axis image;
        axis ij;
        set(handles.activity_axis,'XLim',xlim);
        set(handles.activity_axis,'YLim',ylim);
        axis off;

        if handles.NeedSpeed
            Y=zeros(1,handles.NUM_ELEMENTS-lag);
            for k=1:(handles.NUM_ELEMENTS-lag)
                ind=(handles.Speed{k,lag}>handles.Filt);
                Y(k)=100*sum(ind)/length(ind);
            end 
            tit='Percent of perimeter protruding faster than the Filter value';
        else

            if(get(handles.pro_act,'Value') == 1)
                if(get(handles.perim_norm, 'Value') == 1)            
                    Y = handles.PrActR{lag}(4,:);   
                    tit='Protrusive area normalized by cell perimeter';
                elseif (get(handles.area_norm, 'Value') == 1)    
                    Y = handles.PrActR{lag}(1,:);
                    tit='Protrusive area normalized by cell area';
                elseif (get(handles.none_norm, 'Value') == 1)
                    Y = handles.PrActF{lag};
                    tit='Protrusive area (not normalized)';
                end
            elseif(get(handles.retro_act,'Value') == 1)
                if(get(handles.perim_norm, 'Value') == 1)
                    Y = handles.PrActR{lag}(5,:);      
                    tit='Retractive area normalized by cell perimeter';
                elseif (get(handles.area_norm, 'Value') == 1)    
                    Y = handles.PrActR{lag}(2,:);
                    tit='Retractive area normalized by cell area';
                elseif (get(handles.none_norm, 'Value') == 1)
                    Y = handles.PrActB{lag};
                    tit='Retractive area (not normalized)';
                end  
            elseif(get(handles.total_act,'Value') == 1)
                if(get(handles.perim_norm, 'Value') == 1)                
                    Y = handles.PrActR{lag}(6,:);      
                    tit='Total area normalized by cell perimeter';
                elseif (get(handles.area_norm, 'Value') == 1)    
                    Y = handles.PrActR{lag}(3,:);
                    tit='Total area normalized by cell area';
                elseif (get(handles.none_norm, 'Value') == 1)
                    Y = handles.PrActT{lag};
                    tit='Total area (not normalized)';
                end
            end        
        end

         if get(handles.checkbox2,'Value')
            Z=GFilterA(Y,handles.GN);
        else
            Z=Y;
        end
        [mv,mi]=max(Z);

        set(handles.mean_text,'String', ['Mean = ' num2str(mean(Y),'%.2f')]);
        set(handles.max_text, 'String', ['Max = ' num2str(mv,'%.2f')]);
        set(handles.peak_text,'String', ['at T = ' num2str(mi)]);
        set(handles.result_title,'String',tit);  

        axes(handles.results_axis);
        cla;
        hold on;
        plot([1 handles.NUM_ELEMENTS],[mv mv],':','Color',[.6 .6 .6]);
        plot([mi mi],[0 mv],':','Color',[.6 .6 .6]);
        plot([1 handles.NUM_ELEMENTS],[mean(Y) mean(Y)],'r','LineWidth',2);
        plot(Y,'g','LineWidth',2);
        plot(Z,'b','LineWidth',2);
        if min(Y)<max(Y)
            set(handles.results_axis,'Box','on','XLim',[1 handles.NUM_ELEMENTS],'YLim',[1.1*min(Y)-0.1*max(Y) 1.1*max(Y)-0.1*min(Y)]);
        elseif max(Y)>0
            set(handles.results_axis,'Box','on','XLim',[1 handles.NUM_ELEMENTS],'YLim',[0.9*max(Y) 1.1*max(Y)]);
        else
            set(handles.results_axis,'Box','on','XLim',[1 handles.NUM_ELEMENTS],'YLim',[-1 1]);
        end

        guidata(hObject,handles);
    end   
    set(handle(jWindow),'Enabled',true);
    
function filterSlider_Callback(hObject, eventdata, handles)

    if handles.Done2(handles.Lag)
        
        filt = round(get(handles.filterSlider,'Value'));
        
        Y=zeros(1,handles.NUM_ELEMENTS-handles.Lag);
        for k=1:(handles.NUM_ELEMENTS-handles.Lag)
            ind=(handles.Speed{k,handles.Lag}>filt);
            Y(k)=100*sum(ind)/length(ind);
        end
        
        if get(handles.checkbox2,'Value')
            Z=GFilterA(Y,handles.GN);
        else
            Z=Y;
        end
        [mv,mi]=max(Z);
        set(handles.mean_text,'String', ['Mean = ' num2str(mean(Y),'%.2f')]);
        set(handles.max_text, 'String', ['Max = ' num2str(mv,'%.2f')]);
        set(handles.peak_text,'String', ['at T = ' num2str(mi)]);
        
        axes(handles.results_axis);
        cla;
        hold on;
        plot([1 handles.NUM_ELEMENTS],[mv mv],':','Color',[.6 .6 .6]);
        plot([mi mi],[0 mv],':','Color',[.6 .6 .6]);
        plot([1 handles.NUM_ELEMENTS],[mean(Y) mean(Y)],'r','LineWidth',2);
        plot(Y,'g','LineWidth',2);
        plot(Z,'b','LineWidth',2);
        if min(Y)<max(Y)
            set(handles.results_axis,'Box','on','XLim',[1 handles.NUM_ELEMENTS],'YLim',[1.1*min(Y)-0.1*max(Y) 1.1*max(Y)-0.1*min(Y)]);
        elseif max(Y)>0
            set(handles.results_axis,'Box','on','XLim',[1 handles.NUM_ELEMENTS],'YLim',[0.9*max(Y) 1.1*max(Y)]);
        else
            set(handles.results_axis,'Box','on','XLim',[1 handles.NUM_ELEMENTS],'YLim',[-1 1]);
        end            
        
    end

    
function movieSlider_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end

function lagSlider_CreateFcn(hObject, eventdata, handles)

    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
    
function filterSlider_CreateFcn(hObject, eventdata, handles)

    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
%--------------------------------------------------------------------------

    
  
% --- Buttons -------------------------------------------------------------
function run_lag_range_Callback(hObject, eventdata, handles)
    
    %set(handles.run_lag_range, 'Enable', 'off');
    jFigPeer = get(handle(gcf),'JavaFrame');
    try
        jWindow = jFigPeer.fFigureClient.getWindow;
    catch
        jWindow = jFigPeer.fHG1Client.getWindow;
    end
    set(handle(jWindow),'Enabled',false);
    try   
        if ~handles.NeedSpeed  
            oldlag=handles.Lag;
            for lag=1:handles.Range
                set(handles.lagSlider, 'Value', lag);
                pause(0.001);
                handles.Lag=lag;        
                if handles.Done1(lag)==0
                    handles=LAG_STAT(handles);            
                end
            end
            handles.Lag=oldlag;
            set(handles.lagSlider,'Value',oldlag);
            range=handles.Range;
        else
            range=handles.Range2;
        end

        if handles.NeedSpeed

            Y1=zeros(1,handles.Range2);
            Y2=zeros(1,handles.Range2);
            Z2=zeros(1,handles.Range2);
            TY2=zeros(1,handles.Range2);
            TZ2=zeros(1,handles.Range2);

            oldfilt=handles.Filt;
            for filt=1:handles.Range2
                set(handles.filterSlider, 'Value', filt);
                pause(0.001);
                Y=zeros(1,handles.NUM_ELEMENTS-handles.Lag); 
                for k=1:(handles.NUM_ELEMENTS-handles.Lag)
                    ind=(handles.Speed{k,handles.Lag}>(filt-1));
                    Y(k)=100*sum(ind)/length(ind);
                end
                Z=GFilterA(Y,handles.GN);
                Y1(filt)=mean(Y);
                [Y2(filt),TY2(filt)]=max(Y);       
                [Z2(filt),TZ2(filt)]=max(Z);            
            end
            handles.MeanVelPer=Y1;
            handles.MaxVelPer=Y2;
            handles.zMaxVelPer=Z2;
            handles.TVelPer=TY2;
            handles.zTVelPer=TZ2;            
            
            handles.Filt=oldfilt;
            set(handles.filterSlider,'Value',oldfilt);
            tit='Percent of perimeter protruding faster than the Filter value';
            ylb='percentile (%)';
            xlb='filter (pixel/frame)'; 
            ytick=0:(range-1);
            vw=[-60 10];
        else    
            if(get(handles.pro_act,'Value') == 1)
                if(get(handles.perim_norm, 'Value') == 1)            
                    Y1 = handles.MeanActR(4,:);            
                    Y2 = handles.MaxActR(4,:);
                    Z2 = handles.zMaxActR(4,:);
                    tit='Protrusive area normalized by cell perimeter';
                    ylb='length (pixels)';
                elseif (get(handles.area_norm, 'Value') == 1)    
                    Y1 = handles.MeanActR(1,:);
                    Y2 = handles.MaxActR(1,:);
                    Z2 = handles.zMaxActR(1,:);
                    tit='Protrusive area normalized by cell area';
                    ylb='percentile (%)';
                elseif (get(handles.none_norm, 'Value') == 1)
                    Y1 = handles.MeanActF;
                    Y2 = handles.MaxActF;
                    Z2 = handles.zMaxActF;
                    tit='Protrusive area (not normalized)';
                    ylb='area (pixel^2)';
                end
            elseif(get(handles.retro_act,'Value') == 1)
                if(get(handles.perim_norm, 'Value') == 1)
                    Y1 = handles.MeanActR(5,:);            
                    Y2 = handles.MaxActR(5,:); 
                    Z2 = handles.zMaxActR(5,:); 
                    tit='Retractive area normalized by cell perimeter';
                    ylb='length (pixels)';
                elseif (get(handles.area_norm, 'Value') == 1)    
                    Y1 = handles.MeanActR(2,:);
                    Y2 = handles.MaxActR(2,:);
                    Z2 = handles.zMaxActR(2,:);
                    tit='Retractive area normalized by cell area';
                    ylb='percentile (%)';
                elseif (get(handles.none_norm, 'Value') == 1)
                    Y1 = handles.MeanActB;
                    Y2 = handles.MaxActB;
                    Z2 = handles.zMaxActB;
                    tit='Retractive area (not normalized)';
                    ylb='area (pixel^2)';
                end  
            elseif(get(handles.total_act,'Value') == 1)
                if(get(handles.perim_norm, 'Value') == 1)
                    Y1 = handles.MeanActR(6,:);            
                    Y2 = handles.MaxActR(6,:);          
                    Z2 = handles.zMaxActR(6,:);          
                    tit='Total activity area normalized by cell perimeter';
                    ylb='length (pixels)';
                elseif (get(handles.area_norm, 'Value') == 1)    
                    Y1 = handles.MeanActR(3,:);
                    Y2 = handles.MaxActR(3,:);
                    Z2 = handles.zMaxActR(3,:);
                    tit='Total activity area normalized by cell area';
                    ylb='percentile (%)';
                elseif (get(handles.none_norm, 'Value') == 1)
                    Y1 = handles.MeanActT;
                    Y2 = handles.MaxActT;
                    Z2 = handles.zMaxActT;
                    tit='Total activity area (not normalized)';                
                    ylb='area (pixel^2)';
                end
            end
            xlb='lag (frames)';  
            ytick=1:range;
            vw=[-115 10];
        end

        f=figure; s=subplot(1,1,1); 
        colormap(spring);
        if get(handles.checkbox2,'Value')
            Y=[Y1(1:range)' Z2(1:range)'];
        else
            Y=[Y1(1:range)' Y2(1:range)'];
        end
        bar3(Y,'detached');    
        if max(Y2)>0
            set(s,'Box','off','XLim',[0.5 2.5],'XTick',[],'YLim',[0.5 range+0.5],'YTick',1:range,'YTickLabel',ytick,'ZLim',[0 1.1*max(Y2)],'View',vw);
        else
            set(s,'Box','off','XLim',[0.5 2.5],'XTick',[],'YLim',[0.5 range+0.5],'YTick',1:range,'YTickLabel',ytick,'ZLim',[0 1],'View',vw);
        end
        legend('mean','max','Location','NorthWestOutside');
        title(tit);
        ylabel(xlb); 
        zlabel(ylb);

        %set(handles.run_lag_range, 'Enable', 'on');
        guidata(hObject, handles);
    end
    set(handle(jWindow),'Enabled',true);
    
%--------------------------------------------------------------------------   



% --- Selections ----------------------------------------------------------
function activity_type_SelectionChangeFcn(hObject, eventdata, handles)

    switch get(eventdata.NewValue, 'Tag')
        case 'pro_act'
            if(get(handles.perim_norm, 'Value') == 1)
                Y = handles.PrActR{handles.Lag}(4,:);      
                tit='Protrusive area normalized by cell perimeter';
            elseif (get(handles.area_norm, 'Value') == 1)    
                Y = handles.PrActR{handles.Lag}(1,:);
                tit='Protrusive area normalized by cell area';  
            elseif (get(handles.none_norm, 'Value') == 1)
                Y = handles.PrActF{handles.Lag};
                tit='Protrusive area (not normalized)';
            end
        case 'retro_act'
           if(get(handles.perim_norm, 'Value') == 1)
                Y = handles.PrActR{handles.Lag}(5,:);            
                tit='Retractive area normalized by cell perimeter';
            elseif (get(handles.area_norm, 'Value') == 1)    
                Y = handles.PrActR{handles.Lag}(2,:);
                tit='Retractive area normalized by cell area';
            elseif (get(handles.none_norm, 'Value') == 1)
                Y = handles.PrActB{handles.Lag};
                tit='Retractive area (not normalized)';
            end  
        case 'total_act'
            if(get(handles.perim_norm, 'Value') == 1)
                Y = handles.PrActR{handles.Lag}(6,:); 
                tit='Total area normalized by cell perimeter';
            elseif (get(handles.area_norm, 'Value') == 1)    
                Y = handles.PrActR{handles.Lag}(3,:);
                tit='Total area normalized by cell area';
            elseif (get(handles.none_norm, 'Value') == 1)
                Y = handles.PrActT{handles.Lag};
                tit='Total area (not normalized)';
            end
    end

    set(handles.result_title,'String',tit);
    
    if get(handles.checkbox2,'Value')
        Z=GFilterA(Y,handles.GN);
    else
        Z=Y;
    end
    [mv,mi]=max(Z);
    axes(handles.results_axis);
    cla;
    hold on;
    plot([1 handles.NUM_ELEMENTS],[mv mv],':','Color',[.6 .6 .6]);
    plot([mi mi],[0 mv],':','Color',[.6 .6 .6]);
    plot([1 handles.NUM_ELEMENTS],[mean(Y) mean(Y)],'r','LineWidth',2);
    plot(Y,'g','LineWidth',2);
    plot(Z,'b','LineWidth',2);
    if min(Y)<max(Y)
        set(handles.results_axis,'Box','on','XLim',[1 handles.NUM_ELEMENTS],'YLim',[1.1*min(Y)-0.1*max(Y) 1.1*max(Y)-0.1*min(Y)]);
    elseif max(Y)>0
        set(handles.results_axis,'Box','on','XLim',[1 handles.NUM_ELEMENTS],'YLim',[0.9*max(Y) 1.1*max(Y)]);
    else
        set(handles.results_axis,'Box','on','XLim',[1 handles.NUM_ELEMENTS],'YLim',[-1 1]);
    end
    
    set(handles.mean_text,'String', ['Mean = ' num2str(mean(Y),'%.2f')]);
    set(handles.max_text, 'String', ['Max = ' num2str(mv,'%.2f')]);
    set(handles.peak_text,'String', ['at T = ' num2str(mi)]);

function norm_type_SelectionChangeFcn(hObject, eventdata, handles)
    
     switch get(eventdata.NewValue, 'Tag')
        case 'perim_norm'
            if(get(handles.pro_act, 'Value') == 1)
                Y = handles.PrActR{handles.Lag}(4,:);    
                tit='Protrusive area normalized by cell perimeter';
            elseif (get(handles.retro_act, 'Value') == 1)    
                Y = handles.PrActR{handles.Lag}(5,:);
                tit='Retractive area normalized by cell perimeter';
            elseif (get(handles.total_act, 'Value') == 1)
                Y = handles.PrActR{handles.Lag}(6,:);
                tit='Total area normalized by cell perimeter';
            end
        case 'area_norm'
           if(get(handles.pro_act, 'Value') == 1)
                Y = handles.PrActR{handles.Lag}(1,:);    
                tit='Protrusive area normalized by cell area';  
            elseif (get(handles.retro_act, 'Value') == 1)    
                Y = handles.PrActR{handles.Lag}(2,:);
                tit='Retractive area normalized by cell area';
            elseif (get(handles.total_act, 'Value') == 1)
                Y = handles.PrActR{handles.Lag}(3,:);
                tit='Total area normalized by cell area';
            end  
        case 'none_norm'
            if(get(handles.pro_act, 'Value') == 1)
                Y = handles.PrActF{handles.Lag}; 
                tit='Protrusive area (not normalized)';
            elseif (get(handles.retro_act, 'Value') == 1)    
                Y = handles.PrActB{handles.Lag};
                tit='Retractive area (not normalized)';
            elseif (get(handles.total_act, 'Value') == 1)
                Y = handles.PrActT{handles.Lag};
                tit='Total area (not normalized)';
            end
     end

    set(handles.result_title,'String',tit);
    
    if get(handles.checkbox2,'Value')
        Z=GFilterA(Y,handles.GN);
    else
        Z=Y;
    end
    [mv,mi]=max(Z);
    axes(handles.results_axis);
    cla;
    hold on;    
    plot([1 handles.NUM_ELEMENTS],[mv mv],':','Color',[.6 .6 .6]);
    plot([mi mi],[0 mv],':','Color',[.6 .6 .6]);
    plot([1 handles.NUM_ELEMENTS],[mean(Y) mean(Y)],'r','LineWidth',2);
    plot(Y,'g','LineWidth',2);
    plot(Z,'b','LineWidth',2);
    if min(Y)<max(Y)
        set(handles.results_axis,'Box','on','XLim',[1 handles.NUM_ELEMENTS],'YLim',[1.1*min(Y)-0.1*max(Y) 1.1*max(Y)-0.1*min(Y)]);
    elseif max(Y)>0
        set(handles.results_axis,'Box','on','XLim',[1 handles.NUM_ELEMENTS],'YLim',[0.9*max(Y) 1.1*max(Y)]);
    else
        set(handles.results_axis,'Box','on','XLim',[1 handles.NUM_ELEMENTS],'YLim',[-1 1]);
    end
    
    set(handles.mean_text,'String', ['Mean = ' num2str(mean(Y),'%.2f')]);
    set(handles.max_text, 'String', ['Max = ' num2str(mv,'%.2f')]);
    set(handles.peak_text,'String', ['at T = ' num2str(mi)]);
%--------------------------------------------------------------------------
   


% --- Checkboxes ----------------------------------------------------------
function checkbox1_Callback(hObject, eventdata, handles)

    handles.NeedSpeed=get(hObject,'Value');    
    
    if handles.NeedSpeed
        set(handles.pro_act, 'Enable', 'off');
        set(handles.retro_act, 'Enable', 'off');
        set(handles.total_act, 'Enable', 'off');
        set(handles.perim_norm, 'Enable', 'off');
        set(handles.area_norm, 'Enable', 'off');
        set(handles.none_norm, 'Enable', 'off');
        set(handles.filter_text,'String',['Filter = ' num2str(handles.Filt)],'ForegroundColor','k');
        set(handles.activity_type,'ForegroundColor',[.5 .5 .5]);
        set(handles.norm_type,'ForegroundColor',[.5 .5 .5]);
        set(handles.upto1,'ForegroundColor',[.5 .5 .5]);        
        set(handles.upto2,'ForegroundColor','k');
        set(handles.filterSlider, 'Enable', 'on');
        set(handles.edit_range,'String',' ','Enable', 'off');
        set(handles.edit_range2,'String',num2str(handles.Range2),'Enable', 'on');        
        
    else
        set(handles.pro_act, 'Enable', 'on');
        set(handles.retro_act, 'Enable', 'on');
        set(handles.total_act, 'Enable', 'on');
        set(handles.perim_norm, 'Enable', 'on');
        set(handles.area_norm, 'Enable', 'on');
        set(handles.none_norm, 'Enable', 'on');
        set(handles.filter_text,'String','Filter','ForegroundColor',[.5 .5 .5]);
        set(handles.activity_type,'ForegroundColor','k');
        set(handles.norm_type,'ForegroundColor','k');
        set(handles.upto1,'ForegroundColor','k');
        set(handles.upto2,'ForegroundColor',[.5 .5 .5]); 
        set(handles.filterSlider, 'Enable', 'off');
        set(handles.edit_range,'String',num2str(handles.Range),'Enable', 'on');
        set(handles.edit_range2,'String',' ','Enable', 'off');        
    end
    
    guidata(hObject, handles);
    lagSlider_Callback(hObject, eventdata, handles);
    
function checkbox2_Callback(hObject, eventdata, handles)
    
    if get(hObject,'Value')
        set(handles.edit_smoo,'String',num2str(handles.GN),'Enable','on');
        set(handles.smoo_text,'ForegroundColor','k');
    else
        set(handles.edit_smoo,'String',' ','Enable','off');
        set(handles.smoo_text,'ForegroundColor',[.5 .5 .5]);
    end
    guidata(hObject, handles);
    lagSlider_Callback(hObject, eventdata, handles);
%--------------------------------------------------------------------------



% --- Editing -------------------------------------------------------------
function edit_range_Callback(hObject, eventdata, handles)

    range=round(abs(str2double(get(hObject,'String'))));
    if range<handles.NUM_ELEMENTS
        handles.Range=range;
    end
    set(hObject,'String',num2str(handles.Range));
    guidata(hObject, handles);    
 
function edit_range2_Callback(hObject, eventdata, handles)

    range=round(abs(str2double(get(hObject,'String'))));
    if range<get(handles.filterSlider,'Max')
        handles.Range2=range;
    end
    set(hObject,'String',num2str(handles.Range2));
    guidata(hObject, handles);

function edit_smoo_Callback(hObject, eventdata, handles)

    smoo=round(abs(str2double(get(hObject,'String'))));
    if smoo<handles.NUM_ELEMENTS
        handles.GN=smoo;
    end
    set(hObject,'String',num2str(handles.GN));
    handles.Done1 = zeros(1,handles.NUM_ELEMENTS-1);
    handles.Done2 = zeros(1,handles.NUM_ELEMENTS-1);
    guidata(hObject, handles);
    lagSlider_Callback(hObject, eventdata, handles);
    
function edit_range_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function edit_range2_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end    
    
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

function import_mat_Callback(hObject, eventdata, handles)

    import(hObject,1);

function import_tif_Callback(hObject, eventdata, handles)

    import(hObject,2);
   
function save_as_Callback(hObject, eventdata, handles)

function save_results_Callback(hObject, eventdata, handles)

    Lag=handles.Lag;
    PrActF=handles.PrActF;     PrActB=handles.PrActB;     PrActT=handles.PrActT;     PrActR=handles.PrActR;
    MeanActF=handles.MeanActF; MeanActB=handles.MeanActB; MeanActT=handles.MeanActT; MeanActR=handles.MeanActR;
    MaxActF=handles.MaxActF;   MaxActB=handles.MaxActB;   MaxActT=handles.MaxActT;   MaxActR=handles.MaxActR;
    TActF=handles.TActF;       TActB=handles.TActB;       TActT=handles.TActT;       TActR=handles.TActR;
    zMaxActF=handles.zMaxActF; zMaxActB=handles.zMaxActB; zMaxActT=handles.zMaxActT; zMaxActR=handles.zMaxActR;
    zTActF=handles.zTActF;     zTActB=handles.zTActB;     zTActT=handles.zTActT;     zTActR=handles.zTActR;
    
    MeanVelPer=handles.MeanVelPer; MaxVelPer=handles.MaxVelPer; zMaxVelPer=handles.zMaxVelPer;
    TVelPer=handles.TVelPer;       zTVelPer=handles.zTVelPer;
    
    uisave({'Lag','PrActF','PrActB','PrActT','PrActR',...
           'MeanActF','MeanActB','MeanActT','MeanActR',...
           'MaxActF','MaxActB','MaxActT','MaxActR',...
           'TActF','TActB','TActT','TActR',... 
           'zMaxActF','zMaxActB','zMaxActT','zMaxActR',...
           'zTActF','zTActB','zTActT','zTActR',...
           'MeanVelPer','MaxVelPer','zMaxVelPer','TVelPer','zTVelPer'},'ResultData');
 
function save_images_Callback(hObject, eventdata, handles)
    
    dlg_title = 'Save Images as a TIF-File';
    [FileName,PathName,FilterIndex] = uiputfile('*.tif',dlg_title,'ResultMovie');
    
    if FilterIndex
        
        jFigPeer = get(handle(gcf),'JavaFrame');
        try
            jWindow = jFigPeer.fFigureClient.getWindow;
        catch
            jWindow = jFigPeer.fHG1Client.getWindow;
        end
        set(handle(jWindow),'Enabled',false);
        try        
            imt=(handles.mask{1+handles.Lag}-handles.mask{1}+1)/2;
            cim=zeros([size(imt),3]); 
            cim(:,:,1)=imt; cim(:,:,2)=imt; cim(:,:,3)=imt;
            L=size(imt,1);
            if handles.NeedSpeed
                bound=handles.BOUNDARIES{1+handles.Lag};
                by=bound(:,1)-handles.Y(1)+1; 
                bx=bound(:,2)-handles.X(1)+1;
                ind=(handles.Speed{1,handles.Lag}>handles.Filt);            
                if sum(ind)>0
                    imt(L*(bx(ind)-1)+by(ind))=0;
                    cim(:,:,2)=imt; cim(:,:,3)=imt;                
                end                                      
            end
            imwrite(cim,[PathName FileName],'Compression','none');
            for i=2:(handles.NUM_ELEMENTS-handles.Lag)
                imt=(handles.mask{i+handles.Lag}-handles.mask{i}+1)/2;
                cim(:,:,1)=imt; cim(:,:,2)=imt; cim(:,:,3)=imt;
                if handles.NeedSpeed
                    bound=handles.BOUNDARIES{i+handles.Lag};
                    by=bound(:,1)-handles.Y(1)+1; 
                    bx=bound(:,2)-handles.X(1)+1;
                    ind=(handles.Speed{i,handles.Lag}>handles.Filt);            
                    if sum(ind)>0
                        imt(L*(bx(ind)-1)+by(ind))=0;
                        cim(:,:,2)=imt; cim(:,:,3)=imt;
                    end                                     
                end

                imwrite(cim,[PathName FileName],'Compression','none','WriteMode','append');            
            end
        end    
        set(handle(jWindow),'Enabled',true);
    end
% -------------------------------------------------------------------------
    
    
  
