function [mN,mL,sL]=PlotDetectionResultsX(handles,st)

    ct=handles.Time;
    fw=handles.Width;
    fl=handles.Length;
    mp=handles.MP;

    mL=zeros(1,mp);
    sL=zeros(1,mp);
    mN=zeros(1,mp);
    
    
    axes(handles.progress1);
    cla;
    for fr=1:mp
        bx=[0 fr fr 0 0]/mp;
        by=[0 0 1 1 0];
        fill(bx,by,[.7 .7 .7]);
        set(handles.progress1,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[]);
        pause(0.001);
        
        %{
        Lt=length(handles.Tip{fr});
        filo=cell(1,Lt);
        pd=zeros(1,Lt);

        for i=1:Lt                
            [~,filo{i},pd(i)]=FindPathCond(handles.INDX{fr},handles.Coor{fr},handles.DIST{fr},...
                                           handles.Rad{fr},handles.Vb{fr}(handles.Tip{fr}(i)),fw,handles.CutOff);
        end
        %}
        [~,pd,~]=FindTipPaths(handles.PathInd{fr},handles.PathLng{fr},handles.Coor{fr},handles.DIST{fr},...
                              handles.Rad{fr},handles.Tip{fr},fw,handles.CutOff);
        
        mN(fr)=sum(pd>fl);
        if mN(fr)>0
            mL(fr)=mean(pd(pd>fl));            
            sL(fr)=std(pd(pd>fl));            
        end
    end                              
    cla reset;    
    set(handles.progress1,'Box','on','XTick',[],'YTick',[],'Visible','off');
    
    if st
        
        cla(handles.right_axis);
        set(handles.right_axis,'Visible','off');
        axes(handles.top_axis);
        cla;
        hold on;
        plot(1:mp,mN,'bo-','MarkerFaceColor','b','MarkerSize',4);
        plot(ct,mN(ct),'ro','MarkerFaceColor','r','MarkerSize',6);
        hold off;       
        ylabel('Filopodia Number');
        title(['<N> = ' num2str(mean(mN),'%.2f') ' ' char(177) ' ' num2str(std(mN),'%.2f') ...
               ',  <L> = ' num2str(mean(mL(mL>0)),'%.2f') ' ' char(177) ' ' num2str(std(mL(mL>0)),'%.2f') ...
               '  ( <L> - L_{cr} = ' num2str(mean(mL(mL>0)-fl),'%.2f') ' )']);
        set(handles.top_axis,'Box','on','XLim',[1 mp],'XGrid','on','YGrid','on','Visible','on');%,'FontName','TimesNewRoman');
        
        
        axes(handles.bottom_axis);
        cla;
        hold on;
        x=find(mL>0);        
        plot(x,mL(x),'bo-','MarkerFaceColor','b','MarkerSize',4);
        plot(x,mL(x)-sL(x),'b:');
        plot(x,mL(x)+sL(x),'b:');
        plot([1 mp],[fl fl],'k');
        if mL(ct)>0
            plot(ct,mL(ct),'ro','MarkerFaceColor','r','MarkerSize',6);
        end
        hold off;
        xlabel('Time');
        ylabel('Filopodia Length');
        set(handles.bottom_axis,'Box','on','XLim',[1 mp],'XGrid','on','YGrid','on','Visible','on');
        
    end
