function [POS,LEN,i1,i2,tmtr]=PlotTrackingResultsX(handles,st)

    v1=get(handles.radiobutton1,'Value');
    v2=get(handles.radiobutton2,'Value');

    fw=handles.Width;
    fl=handles.Length;
    fd=handles.Distance;
    fg=handles.Gap;
    ff=handles.Filter;
    mp=handles.MP;

    %cxprs=cell(1,mp);
    %cyprs=cell(1,mp);
    cfilo=cell(1,mp);
    cpd=cell(1,mp);

    POS=zeros(mp,5000);
    LEN=zeros(mp,5000);
    
    axes(handles.progress1);
    cla;
    bx=[0 1 1 0 0]/mp;
    by=[0 0 1 1 0];
    fill(bx,by,[.7 .7 .7]);
    set(handles.progress1,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[]);
    pause(0.001);
        
    %{
    Lt=length(handles.Tip{1});
    filo=cell(1,Lt);
    pd=zeros(1,Lt);    
    for i=1:Lt                
        [~,filo{i},pd(i)]=FindPathCond(handles.INDX{1},handles.Coor{1},handles.DIST{1},...
                                       handles.Rad{1},handles.Vb{1}(handles.Tip{1}(i)),fw,handles.CutOff);
    end
    %}
    [filo,pd,~]=FindTipPaths(handles.PathInd{1},handles.PathLng{1},handles.Coor{1},handles.DIST{1},...
                              handles.Rad{1},handles.Tip{1},fw,handles.CutOff);
    
    filo=filo(pd>fl);
    pd=pd(pd>fl);
    
    cfilo{1}=filo;
    cpd{1}=pd;
    
    P=pd; N=length(P);
    POS(1,1:N)=1:N;
    LEN(1,1:N)=P;
        
    for n=2:mp 
        bx=[0 n n 0 0]/mp;
        by=[0 0 1 1 0];
        fill(bx,by,[.7 .7 .7]);
        set(handles.progress1,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[]);
        pause(0.001);
                 
        %{
        Lt=length(handles.Tip{n});
        filo=cell(1,Lt);
        pd=zeros(1,Lt);
        for i=1:Lt                
            [~,filo{i},pd(i)]=FindPathCond(handles.INDX{n},handles.Coor{n},handles.DIST{n},...
                                           handles.Rad{n},handles.Vb{n}(handles.Tip{n}(i)),fw,handles.CutOff);
        end
        %}
        [filo,pd,~]=FindTipPaths(handles.PathInd{n},handles.PathLng{n},handles.Coor{n},handles.DIST{n},...
                              handles.Rad{n},handles.Tip{n},fw,handles.CutOff);
        filo=filo(pd>fl);
        pd=pd(pd>fl);
        
        cfilo{n}=filo;
        cpd{n}=pd;                                              
        
        filo2=filo;
        P2=pd;
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
            POS(n,(N+1):(N+A))=TIND;
            LEN(n,(N+1):(N+A))=P2(TIND);        
            N=N+A;
        end                

    end
    cla reset;    
    set(handles.progress1,'Box','on','XTick',[],'YTick',[],'Visible','off');
    
    POS=POS(:,1:N);
    LEN=LEN(:,1:N);
    i1=zeros(1,N);
    i2=zeros(1,N);
    tmtr=zeros(1,N);

    if st
        axes(handles.right_axis);
        cla;
        hold on;

        for o=1:N
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

            if o==1            
                line([o o],[i1(o) i2(o)],'Color','y','LineWidth',1);        
            elseif tmtr(o)>ff
                line([o o],[i1(o) i2(o)],'Color','g','LineWidth',1);        
            else
                line([o o],[i1(o) i2(o)],'Color','k','LineWidth',1);        
            end
            if (i2(o)-i1(o))>1
                if o==1 && k
                    plot(o+0*jj,i1(o)+jj,'yo','LineWidth',1,'MarkerFaceColor',[.4 .4 .4],'MarkerSize',4);
                elseif k && tmtr(o)>ff
                    plot(o+0*jj,i1(o)+jj,'go','LineWidth',1,'MarkerFaceColor',[.4 .4 .4],'MarkerSize',4);
                else
                    plot(o+0*jj,i1(o)+jj,'ko','LineWidth',1,'MarkerFaceColor',[.4 .4 .4],'MarkerSize',4);
                end
            end        
            if o==1
                plot([o o],[i1(o) i2(o)],'yo','LineWidth',1,'MarkerFaceColor','y','MarkerSize',4);            
            elseif tmtr(o)>ff
                plot([o o],[i1(o) i2(o)],'go','LineWidth',1,'MarkerFaceColor','g','MarkerSize',4);        
            else
                plot([o o],[i1(o) i2(o)],'ko','LineWidth',1,'MarkerFaceColor','k','MarkerSize',4); 
            end

        end       
        hold off;

        xlabel('Protrusion ID');
        ylabel('Time');
        if ~isempty(i1)
            title(['lasted from  t = ' num2str(i1(1)) '  to t = ' num2str(i2(1)) ' ; tau = ' num2str(tmtr(1))]);
        end        
        set(handles.right_axis,'Box','on','XLim',[0 N+1],'YLim',[0 mp+1],'XGrid','off','YGrid','off','Color',[.4 .4 .4]);
       

    end


