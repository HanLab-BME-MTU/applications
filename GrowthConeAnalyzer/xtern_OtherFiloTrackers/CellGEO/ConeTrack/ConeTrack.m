function varargout = ConeTrack(varargin)
% CONETRACK MATLAB code for ConeTrack.fig
%      CONETRACK, by itself, creates a new CONETRACK or raises the existing
%      singleton*.
%
%      H = CONETRACK returns the handle to a new CONETRACK or the handle to
%      the existing singleton*.
%
%      CONETRACK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONETRACK.M with the given input arguments.
%
%      CONETRACK('Property','Value',...) creates a new CONETRACK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ConeTrack_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ConeTrack_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ConeTrack

% Last Modified by GUIDE v2.5 24-Sep-2013 08:39:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ConeTrack_OpeningFcn, ...
                   'gui_OutputFcn',  @ConeTrack_OutputFcn, ...
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

% --- Executes just before ConeTrack is made visible.
function ConeTrack_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ConeTrack (see VARARGIN)

% Choose default command line output for ConeTrack
    handles.output = hObject;
    def=importdata('Defaults.dat');
    
    handles.Time = 1;
    handles.DefWidth = def(1);
    handles.DefLength = def(2)/2;
    handles.DefDistance = def(3);
    handles.DefGap = def(4);
    handles.DefFilter = def(5);
    handles.Protrusion=1;

    handles.CutOff = 0;
    handles.isfirst = 1;

    handles.shownreslt=0;
    handles.showntrack=0;

    handles.ShGr=0;
    handles.ShPh=0;
    handles.ShTr=0;

    %set(handles.text11,'String',['  ' char(169) ' UNC-CH, Denis Tsygankov'],'ForegroundColor',[.5 .5 .5]);
    set(handles.save,'Enable','off');

    %---------< Detection Settings >--------- 
    set(handles.checkbox1,'Enable','off');
    set(handles.checkbox2,'Enable','off');

    set(handles.slider1,'Enable','off');
    set(handles.slider2,'Enable','off');
    set(handles.slider3,'Enable','off');
    
    set(handles.detection_button,'Enable','off');

    set(handles.text1,'String','disabled','ForegroundColor',[.5 .5 .5]);
    set(handles.text2,'String','        disabled','ForegroundColor',[.5 .5 .5]);
    set(handles.text3,'String','        disabled','ForegroundColor',[.5 .5 .5]);
    set(handles.text4,'String','disabled','ForegroundColor',[.5 .5 .5]);
    set(handles.text5,'String','disabled','ForegroundColor',[.5 .5 .5]);
    set(handles.text5a,'String','disabled','ForegroundColor',[.5 .5 .5]);
    
    set(handles.uipanel1,'ForegroundColor',[.5 .5 .5]);

    %---------< Tracking Settings >--------- 
    set(handles.checkbox3,'Enable','off');

    set(handles.slider4,'Enable','off');
    set(handles.slider5,'Enable','off');
    set(handles.slider6,'Enable','off');
    set(handles.slider7,'Enable','off');

    set(handles.text7,'String','        disabled','ForegroundColor',[.5 .5 .5]);
    set(handles.text8,'String','   disabled','ForegroundColor',[.5 .5 .5]);
    set(handles.text9,'String','    disabled','ForegroundColor',[.5 .5 .5]);
    set(handles.text10,'String','disabled','ForegroundColor',[.5 .5 .5]);

    set(handles.radiobutton1,'Enable','off');
    set(handles.radiobutton2,'Enable','off');
    set(handles.radiobutton3,'Enable','off');

    set(handles.tracking_button,'Enable','off');

    set(handles.uipanel2,'ForegroundColor',[.5 .5 .5]);
    set(handles.uipanel3,'ForegroundColor',[.5 .5 .5]);
    %-----------------------------------------------------------

    handles.sl1=addlistener(handles.slider1,'Value','PostSet',@(src,evnt)mapper1(hObject,src,evnt));
    handles.sl2=addlistener(handles.slider3,'Value','PostSet',@(src,evnt)mapper2(handles,src,evnt));
    handles.sl3=addlistener(handles.slider3,'Value','PostSet',@(src,evnt)mapper3(handles,src,evnt));

    handles.sl4=addlistener(handles.slider4,'Value','PostSet',@(src,evnt)mapper4(handles,src,evnt));
    handles.sl5=addlistener(handles.slider5,'Value','PostSet',@(src,evnt)mapper5(handles,src,evnt));
    handles.sl6=addlistener(handles.slider6,'Value','PostSet',@(src,evnt)mapper6(handles,src,evnt));
    handles.sl7=addlistener(handles.slider7,'Value','PostSet',@(src,evnt)mapper7(handles,src,evnt));

    set(handles.sl2,'Enabled','off');
    set(handles.sl3,'Enabled','off');
    set(handles.sl4,'Enabled','off');
    set(handles.sl5,'Enabled','off');
    set(handles.sl6,'Enabled','off');
    set(handles.sl7,'Enabled','off');

    guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = ConeTrack_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
    try 
        def = [handles.Width; 2*handles.Length; handles.Distance; handles.Gap; handles.Filter];
        id = fopen([cd '/Defaults.dat'],'wt');
        fprintf(id,'%.2f\n',def);
        fclose(id);
    end
    delete(hObject); 
%--------------------------------------------------------------------------



% --- Function for Listeners ----------------------------------------------
function mapper1(hObject,src,evnt) 

    handles = guidata(hObject);
    fr = round(evnt.NewValue);
    
    oldfr = handles.OLDFR; 
    
    if fr~=oldfr

        istracking=get(handles.checkbox3,'Value');
        if ~istracking    
            %fr=round(evnt.NewValue);
            %fr = round(get(handles.slider1,'Value'));          
            
            set(handles.text1,'String',['time = ' num2str(fr)]);

            axes(handles.left_axis);
            XLim=get(handles.left_axis,'XLim');
            YLim=get(handles.left_axis,'YLim');
            cla;
            fill(handles.BOUNDARIES{fr}(:,2),handles.BOUNDARIES{fr}(:,1),[.9 .9 .9]);        
            axis ij;
            axis image;
            axis off;
            set(handles.left_axis,'XLim',XLim,'YLim',YLim);

            if handles.ShPh
                fw=handles.Width;
                fl=handles.Length;
                
                
                %[filo,pd,~]=FindTipPaths(handles.PathInd{fr},handles.PathLng{fr},handles.Coor{fr},handles.DIST{fr},...
                %              handles.Rad{fr},handles.Tip{fr},fw,handles.CutOff);

                PeakInd=handles.sPeakInd{fr};
                Rad=handles.Rad{fr};
                PeakInd=PeakInd(Rad(PeakInd)>fw);
                if fl<fw
                    [ConeBnd,A,cx,cy]=FindConesX(PeakInd,Rad,handles.Coor{fr},handles.PathInd{fr},handles.PathLng{fr},fl);
                    PlotCones(handles,ConeBnd,cx,cy);
                    
                    mN=length(ConeBnd)-1;
                else
                    mN=0;
                end
                if mN>0
                    mL=mean(A(2:end));
                else
                    mL=0;
                end
                
                %mN=sum(pd>fl); 
                %if mN>0
                %    mL=mean(pd(pd>fl));
                %else
                %    mL=0;
                %end
                set(handles.text4,'String',['Total Number = ' num2str(mN,'%d')]);
                set(handles.text5,'String',['Mean Area = ' num2str(round(mL),'%d')]);
                %set(handles.text5a,'String','  ');                

                %PlotPaths(handles,filo,pd,fl);
                

                if handles.shownreslt
                    
                    mp=handles.MP;
                    cla(handles.top_axis);
                    hold(handles.top_axis,'on');
                    plot(handles.top_axis,1:mp,handles.mN,'bo-','MarkerFaceColor','b','MarkerSize',4);
                    plot(handles.top_axis,fr,mN,'ro','MarkerFaceColor','r','MarkerSize',6);
                    ylabel(handles.top_axis,'Cone Number');%,'FontName','TimesNewRoman');
                    title(handles.top_axis,['<N> = ' num2str(mean(handles.mN),'%.2f') ' ' char(177) ' ' num2str(std(handles.mN),'%.2f') ...
                           ',  <L> = ' num2str(mean(handles.mL(handles.mL>0)),'%.2f') ' ' char(177) ' ' num2str(std(handles.mL(handles.mL>0)),'%.2f')]);
                    set(handles.top_axis,'Box','on','XLim',[1 mp],'XGrid','on','YGrid','on','Visible','on');%,'FontName','TimesNewRoman');


                    cla(handles.bottom_axis);
                    hold(handles.bottom_axis,'on');
                    y=handles.mL;
                    z=handles.sL;
                    x=find(y>0);            
                    plot(handles.bottom_axis,x,y(x),'bo-','MarkerFaceColor','b','MarkerSize',4); 
                    plot(handles.bottom_axis,x,y(x)-z(x),'b:'); 
                    plot(handles.bottom_axis,x,y(x)+z(x),'b:'); 
                    plot(handles.bottom_axis,[1 mp],[fl fl],'k');
                    if mL>0
                        plot(handles.bottom_axis,fr,mL,'ro','MarkerFaceColor','r','MarkerSize',6);
                    end
                    xlabel(handles.bottom_axis,'Time');%,'FontName','TimesNewRoman');
                    ylabel(handles.bottom_axis,'Cone Area');%,'FontName','TimesNewRoman');
                    set(handles.bottom_axis,'Box','on','XLim',[1 mp],'XGrid','on','YGrid','on','Visible','on');%,'FontName','TimesNewRoman');        
                end
            end
        else
            fr=round(evnt.NewValue);
            set(handles.text1,'String',['time = ' num2str(fr)]);
           
            fw=handles.Width;
            fl=handles.Length;
            fd=handles.Distance;
            
            axes(handles.left_axis);
            XLim=get(handles.left_axis,'XLim');
            YLim=get(handles.left_axis,'YLim'); 
            cla;
            fill(handles.BOUNDARIES{fr}(:,2),handles.BOUNDARIES{fr}(:,1),[.9 .9 .9],'EdgeColor','k');
            
            PeakInd=handles.sPeakInd{fr};
            Rad=handles.Rad{fr};
            PeakInd=PeakInd(Rad(PeakInd)>fw);
            
            pPeakInd=handles.sPeakInd{fr+1};
            pRad=handles.Rad{fr+1};
            pPeakInd=pPeakInd(pRad(pPeakInd)>fw);
            if fl<fw                
                [ConeBnd,~,cx,cy]=FindConesX(PeakInd,Rad,handles.Coor{fr},handles.PathInd{fr},handles.PathLng{fr},fl);                
                [~,~,pcx,pcy]=FindConesX(pPeakInd,pRad,handles.Coor{fr+1},handles.PathInd{fr+1},handles.PathLng{fr+1},fl);                
                PlotCones(handles,ConeBnd,cx,cy);
                hold on;
                plot(pcx,pcy,'r.');
            end
            
            if handles.showntrack
                fp=handles.Protrusion;
                ConeDist(handles,fr,fw,fl,fd,fp);
            else
                %ConeDist(handles,fr,fw,fl,fd,0);
            end
            axis ij;
            axis image;
            axis off;
            set(handles.left_axis,'XLim',XLim,'YLim',YLim);

        end
    end
    handles.OLDFR=fr; 
    guidata(hObject, handles);

function mapper2(handles,src,evnt) 

    fr=handles.Time;
    fw=round(evnt.NewValue)/10;
    fl=handles.Length;
    
    set(handles.text2,'String',['critical radius = ' num2str(fw,'%.1f')]);    

    axes(handles.left_axis);
    XLim=get(handles.left_axis,'XLim');
    YLim=get(handles.left_axis,'YLim');
    fill(handles.BOUNDARIES{fr}(:,2),handles.BOUNDARIES{fr}(:,1),[.9 .9 .9]);
    axis ij;
    axis image;
    axis off;
    set(handles.left_axis,'XLim',XLim,'YLim',YLim);
    
    PeakInd=handles.sPeakInd{fr};
    Rad=handles.Rad{fr};
    PeakInd=PeakInd(Rad(PeakInd)>fw);
    if fl<fw
        [ConeBnd,A,cx,cy]=FindConesX(PeakInd,Rad,handles.Coor{fr},handles.PathInd{fr},handles.PathLng{fr},fl);
        PlotCones(handles,ConeBnd,cx,cy);
        mN=length(ConeBnd)-1;
    else
        mN=0;
    end
    if mN>0
        mL=mean(A(2:end));
    else
        mL=0;    
    end
    %[filo,pd,~]=FindTipPaths(handles.PathInd{fr},handles.PathLng{fr},handles.Coor{fr},handles.DIST{fr},...
    %                          handles.Rad{fr},handles.Tip{fr},fw,handles.CutOff);

    %mN=sum(pd>fl); 
    %if mN>0
    %    mL=mean(pd(pd>fl));
    %else
    %    mL=0;
    %end
    set(handles.text4,'String',['Total Number = ' num2str(mN,'%d')]);
    set(handles.text5,'String',['Mean Area = ' num2str(round(mL),'%d')]);
    %set(handles.text5a,'String','  ');    
                              
    %PlotPaths(handles,filo,pd,fl);
      
    if handles.shownreslt
        set(handles.detection_button,'String','Show Results');        
        cla(handles.top_axis); title(handles.top_axis,' ')
        cla(handles.bottom_axis);        
    end

function mapper3(handles,src,evnt) 

    fr=handles.Time;
    fw=handles.Width;
    fl=round(evnt.NewValue)/10;

    set(handles.text3,'String',['neck width = ' num2str(2*fl,'%.1f')]);

    axes(handles.left_axis);
    XLim=get(handles.left_axis,'XLim');
    YLim=get(handles.left_axis,'YLim');
    fill(handles.BOUNDARIES{fr}(:,2),handles.BOUNDARIES{fr}(:,1),[.9 .9 .9]);
    axis ij;
    axis image;
    axis off;
    set(handles.left_axis,'XLim',XLim,'YLim',YLim);
    
    PeakInd=handles.sPeakInd{fr};
    Rad=handles.Rad{fr};
    PeakInd=PeakInd(Rad(PeakInd)>fw);
    if fl<fw
        [ConeBnd,A,cx,cy]=FindConesX(PeakInd,Rad,handles.Coor{fr},handles.PathInd{fr},handles.PathLng{fr},fl);
        PlotCones(handles,ConeBnd,cx,cy);
        mN=length(ConeBnd)-1;
    else
        mN=0;
    end
    if mN>0
        mL=mean(A(2:end));
    else
        mL=0;    
    end
    %[filo,pd,~]=FindTipPaths(handles.PathInd{fr},handles.PathLng{fr},handles.Coor{fr},handles.DIST{fr},...
    %                          handles.Rad{fr},handles.Tip{fr},fw,handles.CutOff);

    %mN=sum(pd>fl); 
    %if mN>0
    %    mL=mean(pd(pd>fl));
    %else
    %    mL=0;
    %end
    set(handles.text4,'String',['Total Number = ' num2str(mN,'%d')]);
    set(handles.text5,'String',['Mean Area = ' num2str(round(mL),'%d')]);    

    %PlotPaths(handles,filo,pd,fl);
    
    if handles.shownreslt
        set(handles.detection_button,'String','Show Results');        
        cla(handles.top_axis); title(handles.top_axis,' ');
        cla(handles.bottom_axis);        
    end

function mapper4(handles,src,evnt) 

    fr=handles.Time;
    fw=handles.Width;
    fl=handles.Length;
    fd=round(evnt.NewValue)/10;

    set(handles.text7,'String',['distance = ' num2str(fd,'%.1f')]);

    PeakInd=handles.sPeakInd{fr};
    Rad=handles.Rad{fr};
    PeakInd=PeakInd(Rad(PeakInd)>fw);

    pPeakInd=handles.sPeakInd{fr+1};
    pRad=handles.Rad{fr+1};
    pPeakInd=pPeakInd(pRad(pPeakInd)>fw);
    if fl<fw                
        [ConeBnd,~,cx,cy]=FindConesX(PeakInd,Rad,handles.Coor{fr},handles.PathInd{fr},handles.PathLng{fr},fl);                
        [~,~,pcx,pcy]=FindConesX(pPeakInd,pRad,handles.Coor{fr+1},handles.PathInd{fr+1},handles.PathLng{fr+1},fl);                
    end

    axes(handles.left_axis);
    XLim=get(handles.left_axis,'XLim');
    YLim=get(handles.left_axis,'YLim'); 
    cla;
    fill(handles.BOUNDARIES{fr}(:,2),handles.BOUNDARIES{fr}(:,1),[.9 .9 .9],'EdgeColor','k');

    PlotCones(handles,ConeBnd,cx,cy);
    hold on;
    plot(pcx,pcy,'r.');     
    
    
    %FiloDist(handles,fr,fw,fl,fd,0);
    %ConeDist(handles,fr,fw,fl,fd,0);
    axis ij;
    axis image;
    axis off;
    set(handles.left_axis,'XLim',XLim,'YLim',YLim);
    
    if handles.showntrack 
        cla(handles.right_axis,'reset');
        set(handles.right_axis,'Box','on','XTick',[],'YTick',[]);
        
        set(handles.tracking_button,'String','Show Results');
        set(handles.tracking_button,'Enable','on');
        set(handles.text10,'String','disabled','ForegroundColor',[.5 .5 .5]);
        
        set(handles.save,'Enable','off');
    end

function mapper5(handles,src,evnt) 

    fg=round(evnt.NewValue);
    set(handles.text8,'String',['gap = ' num2str(fg)]);

    if handles.showntrack 
        cla(handles.right_axis,'reset');
        set(handles.right_axis,'Box','on','XTick',[],'YTick',[]);
        
        set(handles.tracking_button,'String','Show Results');
        set(handles.tracking_button,'Enable','on');
        set(handles.text10,'String','disabled','ForegroundColor',[.5 .5 .5]);

        fr=handles.Time;
        fw=handles.Width;
        fl=handles.Length;
        fd=handles.Distance;

        axes(handles.left_axis);   
        %FiloDist(handles,fr,fw,fl,fd,0);   
        %ConeDist(handles,fr,fw,fl,fd,0);
        set(handles.save,'Enable','off');
    
    end

function mapper6(handles,src,evnt) 

    ff=round(evnt.NewValue);
    set(handles.text9,'String',['filter = ' num2str(ff)]);

    if handles.showntrack
        
        tracked=handles.tmtr;
        mtr=mean(tracked(tracked>ff));
        set(handles.text10,'String',['<tau> = ' num2str(mtr,'%.2f')]);
        tmtr=handles.tmtr;
        pr=handles.Protrusion;
        POS=handles.POS;
        i1=handles.i1;
        i2=handles.i2;
        N=size(POS,2);

        axes(handles.right_axis);
        cla;
        hold on;

        for o=1:N
            if o==pr
                line([o o],[i1(o) i2(o)],'Color','y','LineWidth',1);        
            elseif tmtr(o)>ff
                line([o o],[i1(o) i2(o)],'Color','g','LineWidth',1);        
            else    
                line([o o],[i1(o) i2(o)],'Color','k','LineWidth',1);        
            end
            if (i2(o)-i1(o))>1
                ii=(i1(o)+1):(i2(o)-1);
                jj=find(POS(ii,o)==0); k=length(jj);
                if o==pr && k
                    plot(o+0*jj,i1(o)+jj,'yo','LineWidth',1,'MarkerFaceColor',[.4 .4 .4],'MarkerSize',4);
                elseif k && tmtr(o)>ff
                    plot(o+0*jj,i1(o)+jj,'go','LineWidth',1,'MarkerFaceColor',[.4 .4 .4],'MarkerSize',4);
                else
                    plot(o+0*jj,i1(o)+jj,'ko','LineWidth',1,'MarkerFaceColor',[.4 .4 .4],'MarkerSize',4);
                end
            end
            if o==pr && tmtr(o)>ff        
                plot([o o],[i1(o) i2(o)],'yo','LineWidth',1,'MarkerFaceColor','y','MarkerSize',4);        
            elseif o==pr
                plot([o o],[i1(o) i2(o)],'yo','LineWidth',1,'MarkerFaceColor','k','MarkerSize',4);        
            elseif tmtr(o)>ff
                plot([o o],[i1(o) i2(o)],'go','LineWidth',1,'MarkerFaceColor','g','MarkerSize',4);        
            else
                plot([o o],[i1(o) i2(o)],'ko','LineWidth',1,'MarkerFaceColor','k','MarkerSize',4);        
            end

        end

    end

function mapper7(handles,src,evnt) 
    
    fr=handles.Time;
    fw=handles.Width;
    fl=handles.Length;
    fd=handles.Distance;
    ff=handles.Filter;
    pr=round(evnt.NewValue);
    tmtr=handles.tmtr;

    POS=handles.POS;
    i1=handles.i1;
    i2=handles.i2;
    N=size(POS,2);

    axes(handles.right_axis);
    cla;
    hold on;

    for o=1:N
        if o==pr
            line([o o],[i1(o) i2(o)],'Color','y','LineWidth',1);        
        elseif tmtr(o)>ff
            line([o o],[i1(o) i2(o)],'Color','g','LineWidth',1);        
        else    
            line([o o],[i1(o) i2(o)],'Color','k','LineWidth',1);        
        end
        if (i2(o)-i1(o))>1
            ii=(i1(o)+1):(i2(o)-1);
            jj=find(POS(ii,o)==0); k=length(jj);
            if o==pr && k
                plot(o+0*jj,i1(o)+jj,'yo','LineWidth',1,'MarkerFaceColor',[.4 .4 .4],'MarkerSize',4);
            elseif k && tmtr(o)>ff
                plot(o+0*jj,i1(o)+jj,'go','LineWidth',1,'MarkerFaceColor',[.4 .4 .4],'MarkerSize',4);
            else
                plot(o+0*jj,i1(o)+jj,'ko','LineWidth',1,'MarkerFaceColor',[.4 .4 .4],'MarkerSize',4);
            end
        end
        if o==pr && tmtr(o)>ff        
            plot([o o],[i1(o) i2(o)],'yo','LineWidth',1,'MarkerFaceColor','y','MarkerSize',4);        
        elseif o==pr
            plot([o o],[i1(o) i2(o)],'yo','LineWidth',1,'MarkerFaceColor','k','MarkerSize',4);        
        elseif tmtr(o)>ff
            plot([o o],[i1(o) i2(o)],'go','LineWidth',1,'MarkerFaceColor','g','MarkerSize',4);        
        else
            plot([o o],[i1(o) i2(o)],'ko','LineWidth',1,'MarkerFaceColor','k','MarkerSize',4);        
        end

    end

    hold off;

    title(['lasted from  t = ' num2str(i1(pr)) '  to t = ' num2str(i2(pr)) ' ; tua = ' num2str(tmtr(pr))]);

    PeakInd=handles.sPeakInd{fr};
    Rad=handles.Rad{fr};
    PeakInd=PeakInd(Rad(PeakInd)>fw);
    if fl<fw
        [ConeBnd,~,cx,cy]=FindConesX(PeakInd,Rad,handles.Coor{fr},handles.PathInd{fr},handles.PathLng{fr},fl);        
    end
    
    Xb=handles.BOUNDARIES{fr}(:,2);
    Yb=handles.BOUNDARIES{fr}(:,1);
    %pXb=handles.BOUNDARIES{fr+1}(:,2);
    %pYb=handles.BOUNDARIES{fr+1}(:,1);

    axes(handles.left_axis);   
    XLim=get(handles.left_axis,'XLim');
    YLim=get(handles.left_axis,'YLim'); 
    cla;
    fill(Xb,Yb,[.9 .9 .9],'EdgeColor','k');
    PlotCones(handles,ConeBnd,cx,cy);
    hold on;
    ConeDist(handles,fr,fw,fl,fd,pr);
    axis ij;
    axis image;
    axis off;
    set(handles.left_axis,'XLim',XLim,'YLim',YLim);
% -------------------------------------------------------------------------



% --- Additional functions ------------------------------------------------
function import(hObject)
    handles = guidata(hObject);
    
    axes(handles.left_axis);
    colormap(gray); 
    cla; axis fill;    
    set(handles.left_axis,'Box','on','XTick',[],'YTick',[],'Visible','on');
    
    cla(handles.top_axis,'reset');
    set(handles.top_axis,'Visible','off');
    cla(handles.bottom_axis,'reset');
    set(handles.bottom_axis,'Visible','off');
    
    axes(handles.right_axis);
    cla reset; axis fill;    
    set(handles.right_axis,'Box','on','XTick',[],'YTick',[],'Visible','on');
    
    set(handles.save,'Enable','off');    
    
    if ~handles.isfirst
        
        handles.DefWidth = handles.Width;
        handles.DefLength = handles.Length;
        handles.DefDistance = handles.Distance;
        handles.DefGap = handles.Gap;
        handles.DefFilter = handles.Filter;
                
        %---------< Detection Settings >--------- 
        set(handles.checkbox1,'Enable','off');
        set(handles.checkbox2,'Enable','off');

        set(handles.slider1,'Enable','off');
        set(handles.slider2,'Enable','off');
        set(handles.slider3,'Enable','off');
        
        handles.shownreslt = 0;
        set(handles.detection_button,'Enable','off');        
        set(handles.detection_button,'String','Show Results');
        
        set(handles.text1,'String','disabled','ForegroundColor',[.5 .5 .5]);
        set(handles.text2,'String','        disabled','ForegroundColor',[.5 .5 .5]);
        set(handles.text3,'String','        disabled','ForegroundColor',[.5 .5 .5]);
        set(handles.text4,'String','disabled','ForegroundColor',[.5 .5 .5]);
        set(handles.text5,'String','disabled','ForegroundColor',[.5 .5 .5]);
        set(handles.text5a,'String','disabled','ForegroundColor',[.5 .5 .5]);

        set(handles.uipanel1,'ForegroundColor',[.5 .5 .5]);

        %---------< Tracking Settings >--------- 
        set(handles.checkbox3,'Enable','off');

        set(handles.slider4,'Enable','off');
        set(handles.slider5,'Enable','off');
        set(handles.slider6,'Enable','off');
        set(handles.slider7,'Enable','off');

        set(handles.text7,'String','        disabled','ForegroundColor',[.5 .5 .5]);
        set(handles.text8,'String','   disabled','ForegroundColor',[.5 .5 .5]);
        set(handles.text9,'String','    disabled','ForegroundColor',[.5 .5 .5]);
        set(handles.text10,'String','disabled','ForegroundColor',[.5 .5 .5]);

        set(handles.radiobutton1,'Enable','off');
        set(handles.radiobutton2,'Enable','off');
        set(handles.radiobutton3,'Enable','off');

        set(handles.tracking_button,'Enable','off');

        set(handles.uipanel2,'ForegroundColor',[.5 .5 .5]);
        set(handles.uipanel3,'ForegroundColor',[.5 .5 .5]);
        %-----------------------------------------------------------
            
    end    
    
    [filename,pathname,filterindex] = uigetfile('*.mat'); 
    
    if filterindex~=0
        
        sXbYb=[];
        load([pathname filename]);
        handles.filename=filename;
        if ~isempty(sXbYb)
            mp=size(sXbYb,2);       
        else
            mp=0;
        end
        
        if mp>1   
            
            handles.isfirst = 0;
            
            handles.BOUNDARIES=sXbYb; 
            handles.REM1=sREM1;
            handles.REM2=sREM2;
            handles.Nconvex=sNconvex;
            handles.Nconcave=sNconcave;
            handles.Nedge=sNedge;
            
            handles.INDX=sINDX;
            handles.Vb=sVb;
            handles.Coor=sCoor;
            handles.DIST=sDIST;
            handles.Rad=sRad;
            handles.Tip=sTip;
            
            handles.PathInd=sPathInd;
            handles.PathLng=sPathLng;  
            
            if handles.DefWidth > 50
                handles.Width = 50;
            else
                handles.Width = handles.DefWidth;
            end
            if handles.DefLength > 10
                handles.Length = 10;
            else
                handles.Length = handles.DefLength;
            end
                        
            if handles.DefDistance > 50
                handles.Distance = 50;
            else
                handles.Distance = handles.DefDistance;
            end
            
            if handles.DefGap > (mp-1)
                handles.Gap = mp-1;
            else
                handles.Gap = handles.DefGap;
            end
            
            if handles.DefFilter > mp
                handles.Filter = mp;
            else
                handles.Filter = handles.DefFilter;                
            end
            
            set(handles.sl1,'Enabled','off');
            set(handles.sl2,'Enabled','off');
            set(handles.sl3,'Enabled','off');
            set(handles.sl4,'Enabled','off');
            set(handles.sl5,'Enabled','off');
            set(handles.sl6,'Enabled','off');
                        
            fr=handles.Time;
            if fr>mp
                fr=mp;
                handles.Time=mp;
            end
            handles.OLDFR=fr;

            set(handles.slider1,'Min',1);
            set(handles.slider1,'Max',mp);
            set(handles.slider1,'Value',fr);
            set(handles.slider1,'SliderStep',[1/(mp-1) 5/(mp-1)]);

            fw=handles.Width; mw=500;
            set(handles.slider2,'Min',10);
            set(handles.slider2,'Max',mw);
            set(handles.slider2,'Value',10*fw);
            set(handles.slider2,'SliderStep',[1/(mw-5) 5/(mw-5)]);

            fl=handles.Length; ml=100;
            set(handles.slider3,'Min',5);
            set(handles.slider3,'Max',ml);
            set(handles.slider3,'Value',10*fl);
            set(handles.slider3,'SliderStep',[2/(ml-10) 10/(ml-10)]);

            fd=handles.Distance; md=500;
            set(handles.slider4,'Min',1);
            set(handles.slider4,'Max',md);
            set(handles.slider4,'Value',10*fd);
            set(handles.slider4,'SliderStep',[2/(md-1) 10/(md-1)]);

            fg=handles.Gap; mg=mp-1;
            set(handles.slider5,'Min',0);
            set(handles.slider5,'Max',mg);
            set(handles.slider5,'Value',fg);
            set(handles.slider5,'SliderStep',[1/(mg-0) 5/(mg-0)]);

            ff=handles.Filter; mf=mp;
            set(handles.slider6,'Min',0);
            set(handles.slider6,'Max',mf);
            set(handles.slider6,'Value',ff);
            set(handles.slider6,'SliderStep',[1/(mf-0) 5/(mf-0)]);

            handles.MP=mp;           
            
            MIx=zeros(1,mp); MAx=MIx; MIy=MIx; MAy=MIx;            
            
            handles.sPeakInd=cell(1,mp);
            axes(handles.progress1);
            cla;            
            for i=1:mp                
                bx=[0 i i 0 0]/mp;
                by=[0 0 1 1 0];
                fill(bx,by,[.7 .7 .7]);
                set(handles.progress1,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[]);
                pause(0.001);
        
                tmpx=sXbYb{i}(:,2);
                tmpy=sXbYb{i}(:,1);
                INDX=sINDX{i};
                Rad=sRad{i};
                MIx(i)=min(tmpx); MAx(i)=max(tmpx);
                MIy(i)=min(tmpy); MAy(i)=max(tmpy); 
                
                VertNum=length(sRad{i});
                PeakInd=zeros(VertNum,1);
                cnt=0;
                for v=1:VertNum   
                    p1=INDX(INDX(:,1)==v,2);    
                    p2=INDX(INDX(:,2)==v,1);
                    if Rad(v)>max(Rad([p1; p2]))
                        cnt=cnt+1;
                        PeakInd(cnt)=v;                        
                    end
                end
                handles.sPeakInd{i}=PeakInd(1:cnt);                
            end
            cla(handles.progress1);
            set(handles.progress1,'Visible','off');
            handles.xlim=[min(MIx)-10 max(MAx)+10];
            handles.ylim=[min(MIy)-10 max(MAy)+10];
            set(handles.left_axis,'XLim',handles.xlim,'YLim',handles.ylim);
            
            if handles.ShTr
                set(handles.checkbox3,'Value',0);
                set(handles.tracking_button,'String','Show Results');
                handles.showntrack=0;
                handles.ShTr=0;                
            end
            
            if handles.ShTr

                set(handles.slider1,'Min',1);
                set(handles.slider1,'Max',mp-1);
                if fr==mp
                    fr=mp-1; handles.Time=fr;
                end
                set(handles.slider1,'Value',fr);
                set(handles.slider1,'SliderStep',[1/(mp-2) 5/(mp-2)]);
                              

                %---------< Detection Settings >---------     
                set(handles.checkbox1,'Enable','off');
                set(handles.checkbox2,'Enable','off');                
                
                set(handles.slider2,'Enable','off');
                set(handles.slider3,'Enable','off');

                set(handles.text1,'String',['time = ' num2str(fr)]);
                set(handles.text2,'String',['critical radius = ' num2str(fw,'%.1f')],'ForegroundColor','k');
                set(handles.text3,'String',['neck width = ' num2str(2*fl,'%.1f')],'ForegroundColor','k');
                set(handles.text4,'String','disabled','ForegroundColor',[.5 .5 .5]);
                set(handles.text5,'String','disabled','ForegroundColor',[.5 .5 .5]);
                set(handles.text5a,'String','disabled','ForegroundColor',[.5 .5 .5]);

                set(handles.uipanel1,'ForegroundColor',[.5 .5 .5]);                

                %---------< Tracking Settings >--------- 
                
                set(handles.checkbox3,'Enable','on');

                set(handles.slider4,'Enable','on');
                set(handles.slider5,'Enable','on');
                set(handles.slider6,'Enable','on');

                %set(handles.radiobutton1,'Enable','on');
                %set(handles.radiobutton2,'Enable','on');
                %set(handles.radiobutton3,'Enable','on');

                set(handles.tracking_button,'Enable','on');
                set(handles.tracking_button,'String','Show Results');
                handles.showntrack=0;
                
                set(handles.text7,'String',['distance = ' num2str(fd)],'ForegroundColor','k');
                set(handles.text8,'String',['gap = ' num2str(fg)],'ForegroundColor','k');
                set(handles.text9,'String',['filter = ' num2str(ff)],'ForegroundColor','k');

                set(handles.uipanel2,'ForegroundColor','k');
                %set(handles.uipanel3,'ForegroundColor','k');

                %-----------------------------------------------------------
                cla(handles.top_axis,'reset');
                set(handles.top_axis,'Visible','off');
                cla(handles.bottom_axis,'reset');
                set(handles.bottom_axis,'Visible','off');
                cla(handles.right_axis,'reset');
                set(handles.right_axis,'Box','on','XTick',[],'YTick',[],'Visible','on');

                Xb=handles.BOUNDARIES{fr}(:,2);
                Yb=handles.BOUNDARIES{fr}(:,1);
                %pXb=handles.BOUNDARIES{fr+1}(:,2);
                %pYb=handles.BOUNDARIES{fr+1}(:,1);

                axes(handles.left_axis);
                XLim=get(handles.left_axis,'XLim');
                YLim=get(handles.left_axis,'YLim');  
                cla;
                fill(Xb,Yb,[.9 .9 .9],'EdgeColor','k');
                hold on;
                %plot([pXb; pXb(1)],[pYb; pYb(1)],'Color','b');
                %FiloDist(handles,fr,fw,fl,fd,0);        
                axis ij;
                axis image;
                axis off;
                set(handles.left_axis,'XLim',XLim,'YLim',YLim);

                delete(handles.sl4);
                handles.sl4=addlistener(handles.slider4,'Value','PostSet',@(src,evnt)mapper4(handles,src,evnt));

                delete(handles.sl5);
                handles.sl5=addlistener(handles.slider5,'Value','PostSet',@(src,evnt)mapper5(handles,src,evnt));

                delete(handles.sl6);
                handles.sl6=addlistener(handles.slider6,'Value','PostSet',@(src,evnt)mapper6(handles,src,evnt));

            else
            
                set(handles.checkbox1,'Enable','on');
                set(handles.checkbox2,'Enable','on');
                
                set(handles.uipanel1,'ForegroundColor','k');

                axes(handles.left_axis);            
                fill(sXbYb{fr}(:,2),sXbYb{fr}(:,1),[.9 .9 .9]);
                if handles.ShGr
                    PlotGraph(handles,fr);
                end
                axis ij;
                axis image;
                axis off;
                set(handles.left_axis,'XLim',handles.xlim,'YLim',handles.ylim);

                if handles.ShPh

                    %---------< Detection Settings >--------- 
                    set(handles.slider2,'Enable','on');
                    set(handles.slider3,'Enable','on');
                    
                    set(handles.detection_button,'Enable','on');                         
                    
                    set(handles.text2,'String',['critical radius = ' num2str(fw,'%.1f')],'ForegroundColor','k');
                    set(handles.text3,'String',['neck width = ' num2str(2*fl,'%.1f')],'ForegroundColor','k');

                    %---------< Tracking Settings >--------- 
                    set(handles.checkbox3,'Enable','on');

                    set(handles.uipanel2,'ForegroundColor','k');

                    PeakInd=handles.sPeakInd{fr};
                    Rad=handles.Rad{fr};
                    PeakInd=PeakInd(Rad(PeakInd)>fw);
                    if fl<fw
                        [ConeBnd,A,cx,cy]=FindConesX(PeakInd,Rad,handles.Coor{fr},handles.PathInd{fr},handles.PathLng{fr},fl);
                        PlotCones(handles,ConeBnd,cx,cy);
                        mN=length(ConeBnd)-1;
                    else
                        mN=0;
                    end
                    if mN>0
                        mL=mean(A(2:end));
                    else
                        mL=0;    
                    end
                    %{
                    Lt=length(sTip{fr});
                    [filo,pd,tuck]=FindTipPaths(sPathInd{fr},sPathLng{fr},sCoor{fr},sDIST{fr},sRad{fr},sTip{fr},fw,handles.CutOff);
                    
                    nXb=zeros(size(sVb{fr}));
                    nYb=nXb; nLb=0;
                    for i=1:length(sVb{fr})
                        tpath=sPathInd{fr}(1:sPathLng{fr}(i),i);
                        for j=1:Lt
                            tmp=find(tpath==tuck(j),1);
                            if ~isempty(tmp)
                                break;
                            end
                        end
                        if isempty(tmp)
                            nLb=nLb+1;
                            nXb(nLb)=sXbYb{fr}(i,2);
                            nYb(nLb)=sXbYb{fr}(i,1);
                        end
                    end       
                    cXb=[nXb(1:nLb), nXb(1)];
                    cYb=[nYb(1:nLb), nYb(1)];
                    per=sum(sqrt((cXb(2:end)-cXb(1:(end-1))).^2+(cYb(2:end)-cYb(1:(end-1))).^2));                    
                    
                    mN=sum(pd>fl); 
                    if mN>0
                        mL=mean(pd(pd>fl));
                    else
                        mL=0;
                    end
                    %}
                    set(handles.text4,'String',['Total Number = ' num2str(mN,'%d')],'ForegroundColor','k');
                    set(handles.text5,'String',['Mean Area = ' num2str(round(mL),'%d')],'ForegroundColor','k');
                    %set(handles.text5a,'String',['Perimeter = ' num2str(per,'%.2f')],'ForegroundColor','k');
                    
                    %PlotPaths(handles,filo,pd,fl);
                    
                end
                
                delete(handles.sl2);
                handles.sl2=addlistener(handles.slider2,'Value','PostSet',@(src,evnt)mapper2(handles,src,evnt));

                delete(handles.sl3);
                handles.sl3=addlistener(handles.slider3,'Value','PostSet',@(src,evnt)mapper3(handles,src,evnt));
                
            end
            
            set(handles.slider1,'Enable','on');
            set(handles.text1,'String',['time = ' num2str(fr)],'ForegroundColor','k');
                
            delete(handles.sl1);
            handles.sl1=addlistener(handles.slider1,'Value','PostSet',@(src,evnt)mapper1(hObject,src,evnt));
                        
            %set(handles.save,'Enable','on');
        elseif mp==1
            errordlg('The data consists of 1 time frame. The program is meant to analysis movies, not single frames.');
        else
            errordlg('Error is data reading. Please check the file format.');
        end

        
    end
    
    guidata(hObject,handles);       
    
function PlotGraph(handles,time)

    REM1=handles.REM1{time};

    hold(handles.left_axis,'on'); 
    plot(handles.left_axis,[REM1(:,1)'; REM1(:,3)'],[REM1(:,2)';REM1(:,4)'],'Color',[.4 .4 .4]);
    hold(handles.left_axis,'off');

%{    
function PlotPaths(handles,filo,pd,fl)

    hold(handles.left_axis,'on'); 
    for j=1:length(pd)
        if pd(j)>fl
            plot(handles.left_axis,filo{j}(:,1),filo{j}(:,2),'r','LineWidth',2.0);    
        else
            plot(handles.left_axis,filo{j}(:,1),filo{j}(:,2),'y','LineWidth',2.0);    
        end
    end
    hold(handles.left_axis,'off');
%}
    
function PlotCones(handles,ConeBnd,cx,cy)

    hold(handles.left_axis,'on');      
    fill(ConeBnd{1}(:,1),ConeBnd{1}(:,2),[.7 .7 .7]);    
    for j=2:size(ConeBnd,2)
        fill(ConeBnd{j}(:,1),ConeBnd{j}(:,2),'c');    
    end
    plot(cx,cy,'b.');
    hold(handles.left_axis,'off');   
   
function ConeDist(handles,fr,fw,fl,fd,fp)    
    
    cx1=handles.cx{fr};
    cy1=handles.cy{fr};
    cx2=handles.cx{fr+1};
    cy2=handles.cy{fr+1};
    N1=length(cx1);
    N2=length(cx2);
        
    if fp
        N=size(handles.POS,2);
        if N>0
            pos1=handles.POS(fr,fp);    
            ind=find(handles.POS(fr,:)>0);
            tmtr=handles.tmtr;
            ff=handles.Filter;

            for i=ind
                ii=handles.POS(fr,i);
                jj=handles.POS(fr+1,i);

                if ii && jj

                    if ii==pos1 && tmtr(i)>ff
                        plot([cx1(ii) cx2(jj)],[cy1(ii) cy2(jj)],'yo-','MarkerFaceColor','y','MarkerSize',4,'LineWidth',2.0);
                    elseif ii==pos1
                        plot([cx1(ii) cx2(jj)],[cy1(ii) cy2(jj)],'ko-','MarkerFaceColor','k','MarkerSize',4,'LineWidth',2.0);
                    else
                        plot([cx1(ii) cx2(jj)],[cy1(ii) cy2(jj)],'go-','MarkerFaceColor','g','MarkerSize',4,'LineWidth',2.0);            
                    end
                elseif ii 
                    if ii==pos1 && tmtr(i)>ff
                        plot(cx1(ii),cy1(ii),'ro','MarkerFaceColor','y','MarkerSize',4);
                    elseif ii==pos1
                        plot(cx1(ii),cy1(ii),'ko','MarkerFaceColor','y','MarkerSize',4);
                    else
                        plot(cx1(ii),cy1(ii),'ro','MarkerFaceColor','r','MarkerSize',4);
                    end
                end

            end 
        end

    else 
        DM=Inf*ones(N1,N2);

        for i=1:N1         
            for j=1:N2     
                DM(i,j)=(cx1(i)-cx2(j)).*(cx1(i)-cx2(j))+(cy1(i)-cy2(j)).*(cy1(i)-cy2(j));
            end
        end
        BDM=(DM<fd^2);

        for j=1:N2
            if sum(BDM(:,j))>0
                [~,i]=min(DM(:,j));
                 plot([cx1(i) cx2(j)],[cy1(i) cy2(j)],'go-','MarkerFaceColor','g','MarkerSize',4,'LineWidth',2.0);            
            end
        end
    end

%{    
function FiloDist(handles,fr,fw,fl,fd,fp)

    v1=get(handles.radiobutton1,'Value');
    v2=get(handles.radiobutton2,'Value');

    %{
    Lt1=length(handles.Tip{fr});
    filo1=cell(1,Lt1);
    pd1=zeros(1,Lt1);
    for i=1:Lt1                
        [~,filo1{i},pd1(i)]=FindPathCond(handles.INDX{fr},handles.Coor{fr},handles.DIST{fr},...
                                       handles.Rad{fr},handles.Vb{fr}(handles.Tip{fr}(i)),fw,handles.CutOff);
    end
    Lt2=length(handles.Tip{fr+1});
    filo2=cell(1,Lt2);
    pd2=zeros(1,Lt2);

    for i=1:Lt2                
        [~,filo2{i},pd2(i)]=FindPathCond(handles.INDX{fr+1},handles.Coor{fr+1},handles.DIST{fr+1},...
                                       handles.Rad{fr+1},handles.Vb{fr+1}(handles.Tip{fr+1}(i)),fw,handles.CutOff);
    end
    %}
    [filo1,pd1,~]=FindTipPaths(handles.PathInd{fr},handles.PathLng{fr},handles.Coor{fr},handles.DIST{fr},...
                              handles.Rad{fr},handles.Tip{fr},fw,handles.CutOff);
    [filo2,pd2,~]=FindTipPaths(handles.PathInd{fr+1},handles.PathLng{fr+1},handles.Coor{fr+1},handles.DIST{fr+1},...
                              handles.Rad{fr+1},handles.Tip{fr+1},fw,handles.CutOff);
    
    filo1=filo1(pd1>fl);
    filo2=filo2(pd2>fl);
    pd1=pd1(pd1>fl);
    pd2=pd2(pd2>fl);

    N1=length(pd1);
    N2=length(pd2);

    if fp==0
        DM=Inf*ones(N1,N2);

        for i=1:N1         
            for j=1:N2     
                x1=filo1{i}(:,1);
                y1=filo1{i}(:,2);
                n1=length(x1);
                x2=filo2{j}(:,1);
                y2=filo2{j}(:,2);
                n2=length(x2);

                xc1=(x1(1)+x1(end))/2; yc1=(y1(1)+y1(end))/2;
                xc2=(x2(1)+x2(end))/2; yc2=(y2(1)+y2(end))/2;

                dc=abs(xc2-xc1)+abs(yc2-yc1);
                if dc<4*fd                
                    if n1<=n2
                        dis=zeros(1,n1);
                        for k=1:n1
                            dis(k)=min((x1(k)-x2).*(x1(k)-x2)+(y1(k)-y2).*(y1(k)-y2));
                        end
                    else
                        dis=zeros(1,n2);
                        for k=1:n2
                            dis(k)=min((x1-x2(k)).*(x1-x2(k))+(y1-y2(k)).*(y1-y2(k)));
                        end
                    end

                    if v1
                        DM(i,j)=max(dis);
                    elseif v2
                        DM(i,j)=mean(dis);
                    else
                        DM(i,j)=min(dis);
                    end
                end
            end
        end

        BDM=(DM<fd^2);
    end

    for i=1:N1
        plot(filo1{i}(:,1),filo1{i}(:,2),'r','LineWidth',2.0); 
    end    

    for j=1:N2
        plot(filo2{j}(:,1),filo2{j}(:,2),'b','LineWidth',1.0);
    end

    if fp
        N=size(handles.POS,2);
        if N>0
            pos1=handles.POS(fr,fp);    
            ind=find(handles.POS(fr,:)>0);
            tmtr=handles.tmtr;
            ff=handles.Filter;

            for i=ind
                ii=handles.POS(fr,i);
                jj=handles.POS(fr+1,i);

                if ii && jj

                    if ii==pos1 && tmtr(i)>ff
                        plot([filo1{ii}(1,1) filo2{jj}(1,1)],[filo1{ii}(1,2) filo2{jj}(1,2)],'yo-','MarkerFaceColor','y','MarkerSize',4,'LineWidth',2.0);
                    elseif ii==pos1
                        plot([filo1{ii}(1,1) filo2{jj}(1,1)],[filo1{ii}(1,2) filo2{jj}(1,2)],'ko-','MarkerFaceColor','k','MarkerSize',4,'LineWidth',2.0);
                    else
                        plot([filo1{ii}(1,1) filo2{jj}(1,1)],[filo1{ii}(1,2) filo2{jj}(1,2)],'go-','MarkerFaceColor','g','MarkerSize',4,'LineWidth',2.0);            
                    end
                elseif ii 
                    if ii==pos1 && tmtr(i)>ff
                        plot(filo1{ii}(1,1),filo1{ii}(1,2),'ro','MarkerFaceColor','y','MarkerSize',4);
                    elseif ii==pos1
                        plot(filo1{ii}(1,1),filo1{ii}(1,2),'ko','MarkerFaceColor','y','MarkerSize',4);
                    else
                        plot(filo1{ii}(1,1),filo1{ii}(1,2),'ro','MarkerFaceColor','r','MarkerSize',4);
                    end
                end

            end 
        end

    else    

        for j=1:N2
            if sum(BDM(:,j))>0
                [~,i]=min(DM(:,j));
                 plot([filo1{i}(1,1) filo2{j}(1,1)],[filo1{i}(1,2) filo2{j}(1,2)],'go-','MarkerFaceColor','g','MarkerSize',4,'LineWidth',2.0);            
            end
        end
    end
%}    
% -------------------------------------------------------------------------



% --- Sliders -------------------------------------------------------------
function slider1_Callback(hObject, eventdata, handles)

    %fr=round(get(hObject,'Value'));
    fr=round(get(handles.slider1,'Value'));
    handles.Time=fr;
    set(handles.text1,'String',['time = ' num2str(fr)]);
    
    istracking=get(handles.checkbox3,'Value');
    if ~istracking  

        if handles.ShGr
            PlotGraph(handles,fr);
        end

        if handles.ShPh
            fw=handles.Width;
            fl=handles.Length;
            
            PeakInd=handles.sPeakInd{fr};
            Rad=handles.Rad{fr};
            PeakInd=PeakInd(Rad(PeakInd)>fw);
            if fl<fw
                [ConeBnd,~,cx,cy]=FindConesX(PeakInd,Rad,handles.Coor{fr},handles.PathInd{fr},handles.PathLng{fr},fl);
                PlotCones(handles,ConeBnd,cx,cy);
            end
            
            %{
            Lt=length(handles.Tip{fr});            
            [filo,pd,tuck]=FindTipPaths(handles.PathInd{fr},handles.PathLng{fr},handles.Coor{fr},handles.DIST{fr},...
                              handles.Rad{fr},handles.Tip{fr},fw,handles.CutOff);
                          
            nXb=zeros(size(handles.Vb{fr}));
            nYb=nXb; nLb=0;
            for i=1:length(handles.Vb{fr})
                tpath=handles.PathInd{fr}(1:handles.PathLng{fr}(i),i);
                for j=1:Lt
                    tmp=find(tpath==tuck(j),1);
                    if ~isempty(tmp)
                        break;
                    end
                end
                if isempty(tmp)
                    nLb=nLb+1;
                    nXb(nLb)=handles.BOUNDARIES{fr}(i,2);
                    nYb(nLb)=handles.BOUNDARIES{fr}(i,1);
                end
            end       
            cXb=[nXb(1:nLb), nXb(1)];
            cYb=[nYb(1:nLb), nYb(1)];
            per=sum(sqrt((cXb(2:end)-cXb(1:(end-1))).^2+(cYb(2:end)-cYb(1:(end-1))).^2));
            set(handles.text5a,'String',['Perimeter = ' num2str(per,'%.2f')]);
            
            PlotPaths(handles,filo,pd,fl);
            %}
        end    

    end

    delete(handles.sl2);
    handles.sl2=addlistener(handles.slider2,'Value','PostSet',@(src,evnt)mapper2(handles,src,evnt));
    
    delete(handles.sl3);
    handles.sl3=addlistener(handles.slider3,'Value','PostSet',@(src,evnt)mapper3(handles,src,evnt));

    delete(handles.sl4);
    handles.sl4=addlistener(handles.slider4,'Value','PostSet',@(src,evnt)mapper4(handles,src,evnt));

    delete(handles.sl5);
    handles.sl5=addlistener(handles.slider5,'Value','PostSet',@(src,evnt)mapper5(handles,src,evnt));
        
    if handles.showntrack
        delete(handles.sl7);                                
        handles.sl7=addlistener(handles.slider7,'Value','PostSet',@(src,evnt)mapper7(handles,src,evnt));
    end
    
    guidata(hObject, handles);

function slider1_CreateFcn(hObject, eventdata, handles)

    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end

function slider2_Callback(hObject, eventdata, handles)

    fr=handles.Time;
    fw=round(get(hObject,'Value'))/10;
    fl=handles.Length; 

    if handles.ShGr
        PlotGraph(handles,fr);
    end

    if handles.ShPh
        
        PeakInd=handles.sPeakInd{fr};
        Rad=handles.Rad{fr};
        PeakInd=PeakInd(Rad(PeakInd)>fw);
        if fl<fw
            [ConeBnd,~,cx,cy]=FindConesX(PeakInd,Rad,handles.Coor{fr},handles.PathInd{fr},handles.PathLng{fr},fl);
            PlotCones(handles,ConeBnd,cx,cy);
        end

        %{
        Lt=length(handles.Tip{fr});
        [filo,pd,tuck]=FindTipPaths(handles.PathInd{fr},handles.PathLng{fr},handles.Coor{fr},handles.DIST{fr},...
                              handles.Rad{fr},handles.Tip{fr},fw,handles.CutOff);
        nXb=zeros(size(handles.Vb{fr}));
        nYb=nXb; nLb=0;
        for i=1:length(handles.Vb{fr})
            tpath=handles.PathInd{fr}(1:handles.PathLng{fr}(i),i);
            for j=1:Lt
                tmp=find(tpath==tuck(j),1);
                if ~isempty(tmp)
                    break;
                end
            end
            if isempty(tmp)
                nLb=nLb+1;
                nXb(nLb)=handles.BOUNDARIES{fr}(i,2);
                nYb(nLb)=handles.BOUNDARIES{fr}(i,1);
            end
        end       
        cXb=[nXb(1:nLb), nXb(1)];
        cYb=[nYb(1:nLb), nYb(1)];
        per=sum(sqrt((cXb(2:end)-cXb(1:(end-1))).^2+(cYb(2:end)-cYb(1:(end-1))).^2));
        set(handles.text5a,'String',['Perimeter = ' num2str(per,'%.2f')]);
        
        PlotPaths(handles,filo,pd,fl);
        %}
    end

    handles.Width=fw;
    handles.shownreslt=0;

    delete(handles.sl1);
    handles.sl1=addlistener(handles.slider1,'Value','PostSet',@(src,evnt)mapper1(hObject,src,evnt));
    
    delete(handles.sl3);
    handles.sl3=addlistener(handles.slider3,'Value','PostSet',@(src,evnt)mapper3(handles,src,evnt));

    guidata(hObject, handles);
    
function slider2_CreateFcn(hObject, eventdata, handles)

    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end

function slider3_Callback(hObject, eventdata, handles)

    fr=handles.Time;
    fw=handles.Width;
    fl=round(get(hObject,'Value'))/10;

    if handles.ShGr
        PlotGraph(handles,fr);
    end

    if handles.ShPh
        
        PeakInd=handles.sPeakInd{fr};
        Rad=handles.Rad{fr};
        PeakInd=PeakInd(Rad(PeakInd)>fw);
        if fl<fw
            [ConeBnd,~,cx,cy]=FindConesX(PeakInd,Rad,handles.Coor{fr},handles.PathInd{fr},handles.PathLng{fr},fl);
            PlotCones(handles,ConeBnd,cx,cy);
        end
        %[filo,pd,~]=FindTipPaths(handles.PathInd{fr},handles.PathLng{fr},handles.Coor{fr},handles.DIST{fr},...
        %                      handles.Rad{fr},handles.Tip{fr},fw,handles.CutOff);
        %PlotPaths(handles,filo,pd,fl);
    end

    handles.Length=fl;
    handles.shownreslt=0;

    delete(handles.sl1);
    handles.sl1=addlistener(handles.slider1,'Value','PostSet',@(src,evnt)mapper1(hObject,src,evnt));

    delete(handles.sl2);
    handles.sl2=addlistener(handles.slider2,'Value','PostSet',@(src,evnt)mapper2(handles,src,evnt));

    guidata(hObject, handles);

function slider3_CreateFcn(hObject, eventdata, handles)

    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
    
function slider4_Callback(hObject, eventdata, handles)

    fd=round(get(hObject,'Value'))/10;
    set(handles.text7,'String',['distance = ' num2str(fd,'%.1f')]);
    handles.Distance=fd;

    if handles.showntrack 
        handles.showntrack=0;
        set(handles.slider7,'Enable','off');
        set(handles.sl7,'Enabled','off');
    end

    delete(handles.sl1);
    handles.sl1=addlistener(handles.slider1,'Value','PostSet',@(src,evnt)mapper1(hObject,src,evnt));

    delete(handles.sl5);
    handles.sl5=addlistener(handles.slider5,'Value','PostSet',@(src,evnt)mapper5(handles,src,evnt));

    delete(handles.sl6);
    handles.sl6=addlistener(handles.slider6,'Value','PostSet',@(src,evnt)mapper6(handles,src,evnt));

    guidata(hObject, handles);

function slider4_CreateFcn(hObject, eventdata, handles)

    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end

function slider5_Callback(hObject, eventdata, handles)

    fg=round(get(hObject,'Value'));
    set(handles.text8,'String',['gap = ' num2str(fg)]);
    handles.Gap=fg;

    if handles.showntrack 
        
        handles.showntrack=0;
        set(handles.slider7,'Enable','off');
        set(handles.sl7,'Enabled','off');
    end

    delete(handles.sl1);
    handles.sl1=addlistener(handles.slider1,'Value','PostSet',@(src,evnt)mapper1(hObject,src,evnt));

    delete(handles.sl4);
    handles.sl4=addlistener(handles.slider4,'Value','PostSet',@(src,evnt)mapper4(handles,src,evnt));

    delete(handles.sl6);
    handles.sl6=addlistener(handles.slider6,'Value','PostSet',@(src,evnt)mapper6(handles,src,evnt));

    guidata(hObject, handles);

function slider5_CreateFcn(hObject, eventdata, handles)

    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end

function slider6_Callback(hObject, eventdata, handles)

    ff=round(get(hObject,'Value'));
    set(handles.text9,'String',['filter = ' num2str(ff)]);
    handles.Filter=ff;

    delete(handles.sl7);
    handles.sl7=addlistener(handles.slider7,'Value','PostSet',@(src,evnt)mapper7(handles,src,evnt));

    guidata(hObject, handles);

function slider6_CreateFcn(hObject, eventdata, handles)

    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end

function slider7_Callback(hObject, eventdata, handles)

    fp=round(get(hObject,'Value'));
    handles.Protrusion=fp;

    delete(handles.sl1);
    handles.sl1=addlistener(handles.slider1,'Value','PostSet',@(src,evnt)mapper1(hObject,src,evnt));

    delete(handles.sl6);
    handles.sl6=addlistener(handles.slider6,'Value','PostSet',@(src,evnt)mapper6(handles,src,evnt));

    guidata(hObject, handles);

function slider7_CreateFcn(hObject, eventdata, handles)

    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
% -------------------------------------------------------------------------



% --- Buttons -------------------------------------------------------------
function detection_button_Callback(hObject, eventdata, handles)

jFigPeer = get(handle(gcf),'JavaFrame');
try
    jWindow = jFigPeer.fFigureClient.getWindow;
catch
    jWindow = jFigPeer.fHG1Client.getWindow;
end
%set(handle(jWindow),'Enabled',false);

%try
    if handles.shownreslt
        set(handles.detection_button,'String','Show Results');
        handles.shownreslt=0;
        cla(handles.top_axis); title(handles.top_axis,' ')
        cla(handles.bottom_axis);
        
    else        
        [handles.mN,handles.mL,handles.sL]=PlotDetectionResultsY(handles,1);
        set(handles.detection_button,'String','Hide Results');
        handles.shownreslt=1;
    end
      
    delete(handles.sl1);
    handles.sl1=addlistener(handles.slider1,'Value','PostSet',@(src,evnt)mapper1(hObject,src,evnt));
    
    delete(handles.sl2);
    handles.sl2=addlistener(handles.slider2,'Value','PostSet',@(src,evnt)mapper2(handles,src,evnt));
    
    delete(handles.sl3);
    handles.sl3=addlistener(handles.slider3,'Value','PostSet',@(src,evnt)mapper3(handles,src,evnt));
    
    guidata(hObject, handles);    
%end

%set(handle(jWindow),'Enabled',true);
 
function tracking_button_Callback(hObject, eventdata, handles)

jFigPeer = get(handle(gcf),'JavaFrame');
try
    jWindow = jFigPeer.fFigureClient.getWindow;
catch
    jWindow = jFigPeer.fHG1Client.getWindow;
end
%set(handle(jWindow),'Enabled',false);

%try
    fr=handles.Time;
    fw=handles.Width;
    fl=handles.Length;
    fd=handles.Distance;
    ff=handles.Filter;

    if handles.showntrack
        set(handles.tracking_button,'String','Show Results');
        handles.showntrack=0;
        
        cla(handles.right_axis,'reset');
        set(handles.right_axis,'Box','on','XTick',[],'YTick',[]);
        set(handles.slider7,'Enable','off');
        set(handles.sl7,'Enabled','off');
        set(handles.text10,'String','disabled','ForegroundColor',[.5 .5 .5]);

        %axes(handles.left_axis);
        %hold on;
        %FiloDist(handles,fr,fw,fl,fd,0);
        %hold off;
        set(handles.save,'Enable','off');

    else
        set(handles.tracking_button,'String','Hide Results');
        handles.showntrack=1;
        
        %[handles.POS,handles.LEN,handles.i1,handles.i2,handles.tmtr]=PlotTrackingResultsX(handles,1);
        [handles.POS,handles.i1,handles.i2,handles.tmtr,handles.cx,handles.cy]=PlotTrackingResultsY(handles,1);
        
        tracked=handles.tmtr;
        mtr=mean(tracked(tracked>ff));
        set(handles.text10,'String',['<tau> = ' num2str(mtr,'%.2f')],'ForegroundColor','k');
        N=size(handles.POS,2);
        if N>1        
            set(handles.slider7,'Enable','on');
            set(handles.slider7,'Min',1);
            set(handles.slider7,'Max',N);
            set(handles.slider7,'Value',1);
            set(handles.slider7,'SliderStep',[1/(N-1) 5/(N-1)]);    
        end

        handles.Protrusion=1;

        axes(handles.left_axis);
        hold on;
        %FiloDist(handles,fr,fw,fl,fd,1);
        ConeDist(handles,fr,fw,fl,fd,1);
        hold off;
        set(handles.save,'Enable','on');

        delete(handles.sl7);
        handles.sl7=addlistener(handles.slider7,'Value','PostSet',@(src,evnt)mapper7(handles,src,evnt));

    end

    delete(handles.sl1);
    handles.sl1=addlistener(handles.slider1,'Value','PostSet',@(src,evnt)mapper1(hObject,src,evnt));

    delete(handles.sl4);
    handles.sl4=addlistener(handles.slider4,'Value','PostSet',@(src,evnt)mapper4(handles,src,evnt));

    delete(handles.sl5);
    handles.sl5=addlistener(handles.slider5,'Value','PostSet',@(src,evnt)mapper5(handles,src,evnt));

    delete(handles.sl6);
    handles.sl6=addlistener(handles.slider6,'Value','PostSet',@(src,evnt)mapper6(handles,src,evnt));

    guidata(hObject, handles);

%end

%set(handle(jWindow),'Enabled',true);
    
% -------------------------------------------------------------------------    



% --- Checkboxes ----------------------------------------------------------
function checkbox1_Callback(hObject, eventdata, handles)
    
    fr=handles.Time;
    fl=handles.Length;
    handles.ShGr=get(hObject,'Value');
    if handles.ShGr
        PlotGraph(handles,fr);    
    else
        axes(handles.left_axis);
        XLim=get(handles.left_axis,'XLim');
        YLim=get(handles.left_axis,'YLim');
        fill(handles.BOUNDARIES{fr}(:,2),handles.BOUNDARIES{fr}(:,1),[.9 .9 .9]);            
        axis ij;
        axis image;
        axis off;
        set(handles.left_axis,'XLim',XLim,'YLim',YLim);
    end

    if handles.ShPh
        fw=handles.Width;
        
        PeakInd=handles.sPeakInd{fr};
        Rad=handles.Rad{fr};
        PeakInd=PeakInd(Rad(PeakInd)>fw);
        if fl<fw
            [ConeBnd,~,cx,cy]=FindConesX(PeakInd,Rad,handles.Coor{fr},handles.PathInd{fr},handles.PathLng{fr},fl);
            PlotCones(handles,ConeBnd,cx,cy);
        end

        %[filo,pd,~]=FindTipPaths(handles.PathInd{fr},handles.PathLng{fr},handles.Coor{fr},handles.DIST{fr},...
        %                      handles.Rad{fr},handles.Tip{fr},fw,handles.CutOff);
        %PlotPaths(handles,filo,pd,fl);
    end

    guidata(hObject, handles);

function checkbox2_Callback(hObject, eventdata, handles)

    handles.ShPh=get(hObject,'Value');
    fr=handles.Time;
    fw=handles.Width;
    fl=handles.Length;

    if handles.ShPh
        
        PeakInd=handles.sPeakInd{fr};
        Rad=handles.Rad{fr};
        PeakInd=PeakInd(Rad(PeakInd)>fw);
        if fl<fw
            [ConeBnd,A,cx,cy]=FindConesX(PeakInd,Rad,handles.Coor{fr},handles.PathInd{fr},handles.PathLng{fr},fl);
            PlotCones(handles,ConeBnd,cx,cy);
            mN=length(ConeBnd)-1;
        else
            mN=0;
        end
        if mN>0
            mL=mean(A(2:end));
        else
            mL=0;    
        end

    
        %{
        Lt=length(handles.Tip{fr});
        [filo,pd,tuck]=FindTipPaths(handles.PathInd{fr},handles.PathLng{fr},handles.Coor{fr},handles.DIST{fr},...
                              handles.Rad{fr},handles.Tip{fr},fw,handles.CutOff);
        nXb=zeros(size(handles.Vb{fr}));
        nYb=nXb; nLb=0;
        for i=1:length(handles.Vb{fr})
            tpath=handles.PathInd{fr}(1:handles.PathLng{fr}(i),i);
            for j=1:Lt
                tmp=find(tpath==tuck(j),1);
                if ~isempty(tmp)
                    break;
                end
            end
            if isempty(tmp)
                nLb=nLb+1;
                nXb(nLb)=handles.BOUNDARIES{fr}(i,2);
                nYb(nLb)=handles.BOUNDARIES{fr}(i,1);
            end
        end       
        cXb=[nXb(1:nLb), nXb(1)];
        cYb=[nYb(1:nLb), nYb(1)];
        per=sum(sqrt((cXb(2:end)-cXb(1:(end-1))).^2+(cYb(2:end)-cYb(1:(end-1))).^2));

        mN=sum(pd>fl); 
        if mN>0
            mL=mean(pd(pd>fl));
        else
            mL=0;
        end
        PlotPaths(handles,filo,pd,fl);
        %}

        %---------< Detection Settings >---------     
        set(handles.slider2,'Enable','on');
        set(handles.slider3,'Enable','on');
        
        set(handles.detection_button,'Enable','on');

        set(handles.text2,'String',['critical radius = ' num2str(fw,'%.1f')],'ForegroundColor','k');
        set(handles.text3,'String',['neck width = ' num2str(2*fl,'%.1f')],'ForegroundColor','k');
        set(handles.text4,'String',['Total Number = ' num2str(mN,'%d')],'ForegroundColor','k');
        set(handles.text5,'String',['Mean Area = ' num2str(round(mL),'%d')],'ForegroundColor','k');
        %set(handles.text5a,'String',['Perimeter = ' num2str(per,'%.2f')],'ForegroundColor','k');

        %---------< Tracking Settings >--------- 
        set(handles.checkbox3,'Enable','on');

        set(handles.uipanel2,'ForegroundColor','k');

        %-----------------------------------------------------------
        %[handles.mN,handles.mL,handles.sL]=PlotDetectionResultsX(handles,1);
        %[handles.mN,handles.mL,handles.sL]=PlotDetectionResults(handles,1);
        %handles.shownreslt=1;
    else 

        %---------< Detection Settings >---------     
        set(handles.sl2,'Enabled','off');
        set(handles.sl3,'Enabled','off');

        set(handles.slider2,'Enable','off');
        set(handles.slider3,'Enable','off');
        
        handles.shownreslt = 0;
        set(handles.detection_button,'Enable','off');
        set(handles.detection_button,'String','Show Results');

        set(handles.text2,'String','        disabled','ForegroundColor',[.5 .5 .5]);
        set(handles.text3,'String','        disabled','ForegroundColor',[.5 .5 .5]);
        set(handles.text4,'String','disabled','ForegroundColor',[.5 .5 .5]);
        set(handles.text5,'String','disabled','ForegroundColor',[.5 .5 .5]);
        set(handles.text5a,'String','disabled','ForegroundColor',[.5 .5 .5]);

        %---------< Tracking Settings >--------- 
        set(handles.checkbox3,'Enable','off');

        set(handles.uipanel2,'ForegroundColor',[.5 .5 .5]);

        %-----------------------------------------------------------

        axes(handles.left_axis);
        XLim=get(handles.left_axis,'XLim');
        YLim=get(handles.left_axis,'YLim');
        fill(handles.BOUNDARIES{fr}(:,2),handles.BOUNDARIES{fr}(:,1),[.9 .9 .9]);    
        axis ij;
        axis image;
        axis off;
        set(handles.left_axis,'XLim',XLim,'YLim',YLim);

        if handles.ShGr
            PlotGraph(handles,fr);
        end
        
        cla(handles.top_axis,'reset');
        set(handles.top_axis,'Visible','off');
        cla(handles.bottom_axis,'reset');
        set(handles.bottom_axis,'Visible','off');
        set(handles.right_axis,'Visible','on');
        
        %cla(handles.right_axis);
        %handles.needupdate1=0;
        handles.shownreslt=0;
    end

    delete(handles.sl1);
    handles.sl1=addlistener(handles.slider1,'Value','PostSet',@(src,evnt)mapper1(hObject,src,evnt));

    delete(handles.sl2);
    handles.sl2=addlistener(handles.slider2,'Value','PostSet',@(src,evnt)mapper2(handles,src,evnt));

    delete(handles.sl3);
    handles.sl3=addlistener(handles.slider3,'Value','PostSet',@(src,evnt)mapper3(handles,src,evnt));

    guidata(hObject, handles);

function checkbox3_Callback(hObject, eventdata, handles)

    handles.ShTr=get(hObject,'Value');
    fr=handles.Time;
    fw=handles.Width;
    fl=handles.Length;
    fd=handles.Distance;
    fg=handles.Gap;
    ff=handles.Filter;
    mp=handles.MP;

    if handles.ShTr

        set(handles.slider1,'Min',1);
        set(handles.slider1,'Max',mp-1);
        if fr==mp
            fr=mp-1; handles.Time=fr;
        end
        set(handles.slider1,'Value',fr);
        set(handles.slider1,'SliderStep',[1/(mp-2) 5/(mp-2)]);
        
        %---------< Detection Settings >---------     
        set(handles.checkbox1,'Enable','off');
        set(handles.checkbox2,'Enable','off');
       
        set(handles.slider2,'Enable','off');
        set(handles.slider3,'Enable','off');
        handles.shownreslt=0;
        set(handles.detection_button,'Enable','off');
        set(handles.detection_button,'String','Show Results');

        set(handles.text1,'String',['time = ' num2str(fr)]);
        set(handles.text4,'String','disabled','ForegroundColor',[.5 .5 .5]);
        set(handles.text5,'String','disabled','ForegroundColor',[.5 .5 .5]);
        set(handles.text5a,'String','disabled','ForegroundColor',[.5 .5 .5]);

        set(handles.uipanel1,'ForegroundColor',[.5 .5 .5]);

        %---------< Tracking Settings >--------- 

        set(handles.slider4,'Enable','on');
        set(handles.slider5,'Enable','on');
        set(handles.slider6,'Enable','on');

        %set(handles.radiobutton1,'Enable','on');
        %set(handles.radiobutton2,'Enable','on');
        %set(handles.radiobutton3,'Enable','on');
        
        set(handles.tracking_button,'Enable','on');

        set(handles.text7,'String',['distance = ' num2str(fd)],'ForegroundColor','k');
        set(handles.text8,'String',['gap = ' num2str(fg)],'ForegroundColor','k');
        set(handles.text9,'String',['filter = ' num2str(ff)],'ForegroundColor','k');

        %set(handles.uipanel3,'ForegroundColor','k');

        %-----------------------------------------------------------
        cla(handles.top_axis,'reset');
        set(handles.top_axis,'Visible','off');
        cla(handles.bottom_axis,'reset');
        set(handles.bottom_axis,'Visible','off');
        cla(handles.right_axis,'reset');
        set(handles.right_axis,'Box','on','XTick',[],'YTick',[],'Visible','on');

        axes(handles.left_axis);
        XLim=get(handles.left_axis,'XLim');
        YLim=get(handles.left_axis,'YLim'); 
        cla;
        fill(handles.BOUNDARIES{fr}(:,2),handles.BOUNDARIES{fr}(:,1),[.9 .9 .9],'EdgeColor','k');

        PeakInd=handles.sPeakInd{fr};
        Rad=handles.Rad{fr};
        PeakInd=PeakInd(Rad(PeakInd)>fw);

        pPeakInd=handles.sPeakInd{fr+1};
        pRad=handles.Rad{fr+1};
        pPeakInd=pPeakInd(pRad(pPeakInd)>fw);
        if fl<fw                
            [ConeBnd,~,cx,cy]=FindConesX(PeakInd,Rad,handles.Coor{fr},handles.PathInd{fr},handles.PathLng{fr},fl);                
            [~,~,pcx,pcy]=FindConesX(pPeakInd,pRad,handles.Coor{fr+1},handles.PathInd{fr+1},handles.PathLng{fr+1},fl);                
            PlotCones(handles,ConeBnd,cx,cy);
            hold on;
            plot(pcx,pcy,'r.');            
        end
        
        
        %FiloDist(handles,fr,fw,fl,fd,0);  
        
        axis ij;
        axis image;
        axis off;
        set(handles.left_axis,'XLim',XLim,'YLim',YLim);

        delete(handles.sl4);
        handles.sl4=addlistener(handles.slider4,'Value','PostSet',@(src,evnt)mapper4(handles,src,evnt));

        delete(handles.sl5);
        handles.sl5=addlistener(handles.slider5,'Value','PostSet',@(src,evnt)mapper5(handles,src,evnt));

        delete(handles.sl6);
        handles.sl6=addlistener(handles.slider6,'Value','PostSet',@(src,evnt)mapper6(handles,src,evnt));


    else
        
        if handles.showntrack
            cla(handles.right_axis,'reset');
            set(handles.right_axis,'Box','on','XTick',[],'YTick',[]);             
        end

        set(handles.save,'Enable','off');
        
        set(handles.slider1,'Min',1);
        set(handles.slider1,'Max',mp);
        set(handles.slider1,'Value',fr);
        set(handles.slider1,'SliderStep',[1/(mp-1) 5/(mp-1)]);

        %---------< Detection Settings >---------     
        set(handles.checkbox1,'Enable','on');
        set(handles.checkbox2,'Enable','on');

        set(handles.slider1,'Enable','on');
        set(handles.slider2,'Enable','on');
        set(handles.slider3,'Enable','on');
        set(handles.detection_button,'Enable','on');

        set(handles.text1,'String',['time = ' num2str(fr)],'ForegroundColor','k');
        set(handles.text2,'String',['critical radius = ' num2str(fw,'%.1f')],'ForegroundColor','k');
        set(handles.text3,'String',['neck width = ' num2str(2*fl,'%.1f')],'ForegroundColor','k');
       
        set(handles.uipanel1,'ForegroundColor','k');

        %---------< Tracking Settings >--------- 

        set(handles.slider4,'Enable','off');
        set(handles.slider5,'Enable','off');
        set(handles.slider6,'Enable','off');
        set(handles.slider7,'Enable','off');
        set(handles.sl7,'Enabled','off');    

        set(handles.text7,'String','        disabled','ForegroundColor',[.5 .5 .5]);
        set(handles.text8,'String','   disabled','ForegroundColor',[.5 .5 .5]);
        set(handles.text9,'String','    disabled','ForegroundColor',[.5 .5 .5]);
        set(handles.text10,'String','disabled','ForegroundColor',[.5 .5 .5]);

        set(handles.radiobutton1,'Enable','off');
        set(handles.radiobutton2,'Enable','off');
        set(handles.radiobutton3,'Enable','off');

        set(handles.tracking_button,'String','Show Results');
        set(handles.tracking_button,'Enable','off');    
        handles.showntrack=0;

        set(handles.uipanel3,'ForegroundColor',[.5 .5 .5]);

        %-----------------------------------------------------------

        axes(handles.left_axis);
        XLim=get(handles.left_axis,'XLim');
        YLim=get(handles.left_axis,'YLim');
        cla;
        fill(handles.BOUNDARIES{fr}(:,2),handles.BOUNDARIES{fr}(:,1),[.9 .9 .9]);    
        axis ij;
        axis image;
        axis off;
        set(handles.left_axis,'XLim',XLim,'YLim',YLim);
        
        if handles.ShPh
           
            PeakInd=handles.sPeakInd{fr};
            Rad=handles.Rad{fr};
            PeakInd=PeakInd(Rad(PeakInd)>fw);
            if fl<fw
                [ConeBnd,A,cx,cy]=FindConesX(PeakInd,Rad,handles.Coor{fr},handles.PathInd{fr},handles.PathLng{fr},fl);
                PlotCones(handles,ConeBnd,cx,cy);
                mN=length(ConeBnd)-1;                
            else
                mN=0;
            end
            if mN>0
                mL=mean(A(2:end));
            else
                mL=0;    
            end
                
            %{
            Lt=length(handles.Tip{fr});
            [filo,pd,tuck]=FindTipPaths(handles.PathInd{fr},handles.PathLng{fr},handles.Coor{fr},handles.DIST{fr},...
                              handles.Rad{fr},handles.Tip{fr},fw,handles.CutOff);
            nXb=zeros(size(handles.Vb{fr}));
            nYb=nXb; nLb=0;
            for i=1:length(handles.Vb{fr})
                tpath=handles.PathInd{fr}(1:handles.PathLng{fr}(i),i);
                for j=1:Lt
                    tmp=find(tpath==tuck(j),1);
                    if ~isempty(tmp)
                        break;
                    end
                end
                if isempty(tmp)
                    nLb=nLb+1;
                    nXb(nLb)=handles.BOUNDARIES{fr}(i,2);
                    nYb(nLb)=handles.BOUNDARIES{fr}(i,1);
                end
            end       
            cXb=[nXb(1:nLb), nXb(1)];
            cYb=[nYb(1:nLb), nYb(1)];
            per=sum(sqrt((cXb(2:end)-cXb(1:(end-1))).^2+(cYb(2:end)-cYb(1:(end-1))).^2));

            PlotPaths(handles,filo,pd,fl);

            mN=sum(pd>fl); 
            if mN>0
                mL=mean(pd(pd>fl));
            else
                mL=0;
            end
            %}
            set(handles.text4,'String',['Total Number = ' num2str(mN,'%d')],'ForegroundColor','k');
            set(handles.text5,'String',['Mean Area = ' num2str(round(mL),'%d')],'ForegroundColor','k');
            %set(handles.text5a,'String',['Perimeter = ' num2str(per,'%.2f')],'ForegroundColor','k');

            
        end
    end

    delete(handles.sl1);
    handles.sl1=addlistener(handles.slider1,'Value','PostSet',@(src,evnt)mapper1(hObject,src,evnt));
   
    guidata(hObject, handles);
% -------------------------------------------------------------------------



% --- Radiobuttons --------------------------------------------------------
%{
function uipanel3_SelectionChangeFcn(hObject, eventdata, handles)

    fr=handles.Time;
    fw=handles.Width;
    fl=handles.Length;
    fd=handles.Distance;
    
    Xb=handles.BOUNDARIES{fr}(:,2);
    Yb=handles.BOUNDARIES{fr}(:,1);
    pXb=handles.BOUNDARIES{fr+1}(:,2);
    pYb=handles.BOUNDARIES{fr+1}(:,1);

    axes(handles.left_axis);   
    XLim=get(handles.left_axis,'XLim');
    YLim=get(handles.left_axis,'YLim'); 
    cla;
    fill(Xb,Yb,[.9 .9 .9],'EdgeColor','r');
    hold on;
    plot([pXb; pXb(1)],[pYb; pYb(1)],'Color','b');
    FiloDist(handles,fr,fw,fl,fd,0);
    axis ij;
    axis image;
    axis off;
    set(handles.left_axis,'XLim',XLim,'YLim',YLim);

    if handles.showntrack 
        cla(handles.right_axis,'reset');
        set(handles.right_axis,'Box','on','XTick',[],'YTick',[]);
        
        set(handles.tracking_button,'String','Show Results');
        handles.showntrack=0;
    
        set(handles.slider7,'Enable','off');
        set(handles.sl7,'Enabled','off');
        set(handles.text10,'String','disabled','ForegroundColor',[.5 .5 .5]);
    end

    delete(handles.sl1);
    handles.sl1=addlistener(handles.slider1,'Value','PostSet',@(src,evnt)mapper1(hObject,src,evnt));

    delete(handles.sl4);
    handles.sl4=addlistener(handles.slider4,'Value','PostSet',@(src,evnt)mapper4(handles,src,evnt));

    delete(handles.sl5);
    handles.sl5=addlistener(handles.slider5,'Value','PostSet',@(src,evnt)mapper5(handles,src,evnt));

    delete(handles.sl6);
    handles.sl6=addlistener(handles.slider6,'Value','PostSet',@(src,evnt)mapper6(handles,src,evnt));

    guidata(hObject, handles);
%}
% -------------------------------------------------------------------------



% --- Menu ----------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)

function import_Callback(hObject, eventdata, handles)

    import(hObject);
    
function save_Callback(hObject, eventdata, handles)

jFigPeer = get(handle(gcf),'JavaFrame');
try
    jWindow = jFigPeer.fFigureClient.getWindow;
catch
    jWindow = jFigPeer.fHG1Client.getWindow;
end
%set(handle(jWindow),'Enabled',false);

%try
    fw=handles.Width;
    fl=handles.Length;
    fd=handles.Distance;
    fg=handles.Gap;
    ff=handles.Filter;
    mp=handles.MP;
    
    if handles.showntrack
        fp=handles.Protrusion;
        FL=2*fl;
        
        sfilename=['Cone_' num2str(fp,'%02d') '_'  handles.filename];
        
        i1=handles.i1(fp);
        i2=handles.i2(fp);
        
        nConeBnd=cell(1,i2-i1+1);
        nA=zeros(i2-i1+1,1);
        nCx=zeros(i2-i1+1,1);
        nCy=zeros(i2-i1+1,1);
        
        axes(handles.progress1);
        cla;
        for fr=i1:i2
            bx=[0 (fr-i1+1) (fr-i1+1) 0 0]/(i2-i1+1);
            by=[0 0 1 1 0];
            fill(bx,by,[.7 .7 .7]);
            set(handles.progress1,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[]);
            pause(0.001);
        
            n=handles.POS(fr,fp);
            
            PeakInd=handles.sPeakInd{fr};
            Rad=handles.Rad{fr};
            PeakInd=PeakInd(Rad(PeakInd)>fw);
            if fl<fw
                [ConeBnd,A,cx,cy]=FindConesX(PeakInd,Rad,handles.Coor{fr},handles.PathInd{fr},handles.PathLng{fr},fl);                
                nConeBnd{fr-i1+1}=ConeBnd{n};
                nA(fr-i1+1)=A(n);
                nCx(fr-i1+1)=cx(n);
                nCy(fr-i1+1)=cy(n);
            end
                        
        end
        cla(handles.progress1);
        set(handles.progress1,'Visible','off');
        uisave({'fw','FL','fd','fg','ff','mp','fp','i1','i2','nConeBnd','nA','nCx','nCy'},sfilename);  
    end
    
    %{
    N=zeros(1,mp);
    mL=zeros(1,mp);
    sL=zeros(1,mp);
    Per=zeros(1,mp);
    
    axes(handles.progress1);
    cla;
    for fr=1:mp
        bx=[0 fr fr 0 0]/mp;
        by=[0 0 1 1 0];
        fill(bx,by,[.7 .7 .7]);
        set(handles.progress1,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[]);
        pause(0.001);
       
        Lt=length(handles.Tip{fr});
        [~,pd,tuck]=FindTipPaths(handles.PathInd{fr},handles.PathLng{fr},handles.Coor{fr},handles.DIST{fr},...
                              handles.Rad{fr},handles.Tip{fr},fw,handles.CutOff);
        nXb=zeros(size(handles.Vb{fr}));
        nYb=nXb; nLb=0;
        for i=1:length(handles.Vb{fr})
            tpath=handles.PathInd{fr}(1:handles.PathLng{fr}(i),i);
            for j=1:Lt
                tmp=find(tpath==tuck(j),1);
                if ~isempty(tmp)
                    break;
                end
            end
            if isempty(tmp)
                nLb=nLb+1;
                nXb(nLb)=handles.BOUNDARIES{fr}(i,2);
                nYb(nLb)=handles.BOUNDARIES{fr}(i,1);
            end
        end       
        cXb=[nXb(1:nLb), nXb(1)];
        cYb=[nYb(1:nLb), nYb(1)];
        Per(fr)=sum(sqrt((cXb(2:end)-cXb(1:(end-1))).^2+(cYb(2:end)-cYb(1:(end-1))).^2));
        
        N(fr)=sum(pd>fl); 
        if N(fr)>0
            mL(fr)=mean(pd(pd>fl));
            sL(fr)=std(pd(pd>fl));        
        end
    end
    
    cla;
    bx=[0 1 1 0 0]/mp;
    by=[0 0 1 1 0];
    fill(bx,by,[.7 .7 .7]);
    set(handles.progress1,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[]);
    pause(0.001);
    
    v1=get(handles.radiobutton1,'Value');
    v2=get(handles.radiobutton2,'Value');
    %---------------------------------
    cfilo=cell(1,mp);
    cpd=cell(1,mp);

    POS=zeros(mp,5000);
    LEN=zeros(mp,5000);
    
    [filo,pd,~]=FindTipPaths(handles.PathInd{1},handles.PathLng{1},handles.Coor{1},handles.DIST{1},...
                              handles.Rad{1},handles.Tip{1},fw,handles.CutOff);
    filo=filo(pd>fl);
    pd=pd(pd>fl);
    
    cfilo{1}=filo;
    cpd{1}=pd;
    
    P=pd; Npr=length(P);
    POS(1,1:Npr)=1:Npr;
    LEN(1,1:Npr)=P;

    for n=2:mp
        bx=[0 n n 0 0]/mp;
        by=[0 0 1 1 0];
        fill(bx,by,[.7 .7 .7]);
        set(handles.progress1,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[]);
        pause(0.001);
        
        [filo,pd,~]=FindTipPaths(handles.PathInd{n},handles.PathLng{n},handles.Coor{n},handles.DIST{n},...
                              handles.Rad{n},handles.Tip{n},fw,handles.CutOff);
        filo=filo(pd>fl);
        pd=pd(pd>fl);
        
        cfilo{n}=filo;
        cpd{n}=pd;                                              
        
        filo2=filo;
        P2=pd;cx
        N2=length(P2);     

        TIND=1:N2;
        
        for m=(n-1):(-1):max(1,n-fg-1)
            
            filo1=cfilo{m};
            P1=cpd{m};
            N1=length(P1);     

            DM=Inf*ones(N1,N2);

            for i=1:N1

                for j=1:N2            
                    
                    x1=filo1{i}(:,1);
                    y1=filo1{i}(:,2);                    
                    n1=length(x1);
                    x2=filo2{j}(:,1);
                    y2=filo2{j}(:,2);
                    n2=length(x2);

                    xc1=(x1(1)+x1(end))/2; yc1=(y1(1)+y1(end))/2;
                    xc2=(x2(1)+x2(end))/2; yc2=(y2(1)+y2(end))/2;

                    dc=abs(xc2-xc1)+abs(yc2-yc1);
                    if dc<4*fd                
                        if n1<=n2
                            dis=zeros(1,n1);
                            for k=1:n1
                                dis(k)=min((x1(k)-x2).*(x1(k)-x2)+(y1(k)-y2).*(y1(k)-y2));
                            end
                        else
                            dis=zeros(1,n2);
                            for k=1:n2
                                dis(k)=min((x1-x2(k)).*(x1-x2(k))+(y1-y2(k)).*(y1-y2(k)));
                            end
                        end

                        if v1
                            DM(i,j)=max(dis);
                        elseif v2
                            DM(i,j)=mean(dis);
                        else
                            DM(i,j)=min(dis);
                        end
                    end

                end
            end

            BDM=(DM<fd^2);

            temp=ones(1,N2);
            for j=1:N2

                if sum(BDM(:,j))>0
                    [~,i]=min(DM(:,j)); 
                    rind=find(POS(m,:)==i);
                    if POS(n,rind)==0
                        POS(n,rind)=TIND(j);
                        LEN(n,rind)=P2(TIND(j));
                        temp(j)=0;
                    end

                end

            end

            TIND=TIND(temp==1);
            filo2=filo2(temp==1);            
            N2=sum(temp);
            if N2==0
                break;
            end 

        end

        A=sum(temp);
        if A>0
            POS(n,(Npr+1):(Npr+A))=TIND;
            LEN(n,(Npr+1):(Npr+A))=P2(TIND);        
            Npr=Npr+A;
        end                

    end

    POS=POS(:,1:Npr);
    LEN=LEN(:,1:Npr);
    i1=zeros(1,Npr);
    i2=zeros(1,Npr);
    tmtr=zeros(1,Npr);
    tmlg=zeros(1,Npr);
    mxlg=zeros(1,Npr);

    for o=1:Npr
        i1(o)=find(POS(:,o)>0,1,'first');
        i2(o)=find(POS(:,o)>0,1,'last');
        if (i2(o)-i1(o))>1
            ii=(i1(o)+1):(i2(o)-1);
            jj=find(POS(ii,o)==0); k=length(jj);
            tmtr(o)=i2(o)-i1(o)+1-k;
        elseif (i2(o)-i1(o))==1
            tmtr(o)=2;
        else
            tmtr(o)=1;
        end
        if ~isempty(POS(:,o)>0)
            tmlg(o)=mean(LEN(POS(:,o)>0,o));
            mxlg(o)=max(LEN(POS(:,o)>0,o));
        end
    end
    
    mtPer=mean(Per);
    stPer=std(Per);
    mtN=mean(N);
    stN=std(N);
    mtL=mean(mL(mL>0));
    stL=std(mL(mL>0));
    
    prN=sum(tmtr>ff);    
    prT=tmtr(tmtr>ff);
    prL=tmlg(tmtr>ff);
    prML=mxlg(tmtr>ff);
    
    mprT=mean(prT);
    sprT=std(prT);
    mprL=mean(prL);
    sprL=std(prL);
    mprML=mean(prML);
    sprML=std(prML);
    
    crW=fw; crL=2*fl; crD=fd; crG=fg; crF=ff;
    
    cla reset;    
    set(handles.progress1,'Box','on','XTick',[],'YTick',[],'Visible','off');

    uisave({'Per','mtPer','stPer','N','mtN','stN','mL','sL','mtL','stL',...
            'prN','prT','mprT','sprT','prL','mprL','sprL','prML','mprML','sprML',...
            'crW','crL','crD','crG','crF'});  
    %}
%end
%set(handle(jWindow),'Enabled',true);        
        
% -------------------------------------------------------------------------
