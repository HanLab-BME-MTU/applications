function [REM1,REM2,Status,Nconvex,Nconcave,rem]=FindGraphX(Xb,Yb,Lb,handles)

    XB=Xb; YB=Yb; LB=Lb; Status=0;
    MB=[max(XB) min(XB) max(YB) min(YB)];

    AL=zeros(1,LB); AR=zeros(1,LB);
    XO=[]; YO=[]; AO=[]; TO=[]; PO=[];

    Nconvex=0; Nconcave=0;  rem=0; rem1=0; rem2=0;
    for i=1:LB    
        if i==LB
            x0=XB(i-1); x1=XB(i); x2=XB(1);
            y0=YB(i-1); y1=YB(i); y2=YB(1);
        elseif i==1;     
            x0=XB(LB); x1=XB(i); x2=XB(i+1);
            y0=YB(LB); y1=YB(i); y2=YB(i+1);
        else        
            x0=XB(i-1); x1=XB(i); x2=XB(i+1);
            y0=YB(i-1); y1=YB(i); y2=YB(i+1);
        end    

        AR(i)=atan2(y2-y1,x2-x1);
        if AR(i)<0
           AR(i)=AR(i)+2*pi;
        end

        AL(i)=atan2(y0-y1,x0-x1);
        if AL(i)<0
           AL(i)=AL(i)+2*pi;
        end

        pntyp = sign((x2-x1)*(y0-y1) - (x0-x1)*(y2-y1));
        if pntyp==0
            REM1=[]; REM2=[];
            return;
        elseif pntyp==1;
            gam=atan2(sin(AL(i))+sin(AR(i)),cos(AL(i))+cos(AR(i)));
            if gam<0
                gam=gam+2*pi;
            end
            XO=[XO x1];  YO=[YO y1];
            AO=[AO gam]; TO=[TO pntyp];        
            PO=[PO; [x1 y1 0 0 i i]]; 
            Nconvex=Nconvex+1;
        elseif pntyp==-1        
            gam1=AL(i)-pi/2;
            if gam1<0
                gam1=gam1+2*pi;
            end
            gam2=AR(i)+pi/2;
            if gam2>(2*pi)
                gam2=gam2-2*pi;
            end
            XO=[XO x1 x1];     YO=[YO y1 y1];
            AO=[AO gam1 gam2]; TO=[TO pntyp pntyp];
            PO=[PO; [x1 y1 0 0 i i]; [x1 y1 0 0 i i]]; 
            Nconcave=Nconcave+1;
        end   

    end
    LB=length(XO);
    NO=1:LB;
    MEM=NO;
    %--------------------------------------------------------------------------

    %col={'g','r','b','r','b','r','b','r','b','r','b','r','b','r','b','r','b','r','b'};
    
    %if ColCho
    %    col={'b','c'};
    %else
    %    col={'r','m'};
    %end
    
    %--------------------------------------------------------------------------
    Edge=2*(Nconvex+Nconcave)+2*Nconcave-3;
    
    REM1=zeros(2*Nconvex+2*Nconcave-3,4);
    REM2=zeros(2*Nconcave,4);
        
    e=0;
    while LB>2    
        e=e+1;
        
        axes(handles.progress1);
        cla reset;    
        set(handles.progress1,'Box','on','XTick',[],'YTick',[],'Visible','on');
        bx=[0 rem rem 0 0]/Edge;
        by=[0 0 1 1 0];
        fill(bx,by,[.7 .7 .7]);
        set(handles.progress1,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[]);
        pause(0.001);
        
        %axes(handles.boundary_axis);
        %hold on;
        
        [NumInter,Xinter,Yinter,Ainter,Tinter,Pinter]=CollisionA(XB,YB,MB,AL,AR,XO,YO,AO,TO,PO,MEM);
        if isempty(NumInter)             
            MEM=1:length(XO);
            [NumInter,Xinter,Yinter,Ainter,Tinter,Pinter]=CollisionA(XB,YB,MB,AL,AR,XO,YO,AO,TO,PO,MEM);
            if NumInter==0
                %display('NumInter is zero');
                REM1=[]; REM2=[];
                break;
            end
        elseif  NumInter==0 
            MEM=1:length(XO);
            [NumInter,Xinter,Yinter,Ainter,Tinter,Pinter]=CollisionA(XB,YB,MB,AL,AR,XO,YO,AO,TO,PO,MEM);
            if NumInter==0 
                %display('NumInter is zero');
                REM1=[]; REM2=[];
                break;
            end
        elseif NumInter==-1
            REM1=[]; REM2=[];
            break;
        end
        
        LD=length(NumInter);  
        
        if LB==3 && LD>1             
            if max(abs(Xinter-mean(Xinter)))>1e-9 || max(abs(Yinter-mean(Yinter)))>1e-9 
                %display('LB=3, LD>1: Not Fixed');
                REM1=[]; REM2=[];
                break;
            end
            %display('LB=3, LD>1: Fixed');
            LD=1;
            NumInter=NumInter(1);
            Xinter=Xinter(1);
            Yinter=Yinter(1);
            Ainter=Ainter(1);
            Tinter=Tinter(1);
            Pinter=Pinter(1);
        elseif LB==3 && LD==0            
            [NumInter,Xinter,Yinter,Ainter,Tinter,Pinter]=Collision3FIX(XB,YB,MB,AL,AR,XO,YO,AO,TO,PO,MEM);
            
            if length(Xinter)>1 && (max(Xinter)-min(Xinter))<1e-6 && (max(Yinter)-min(Yinter))<1e-6 
                %display('LB=3, LD=0: Fixed');
                LD=1;
                NumInter=NumInter(1);
                Xinter=Xinter(1);
                Yinter=Yinter(1);
                Ainter=Ainter(1);
                Tinter=Tinter(1);
                Pinter=Pinter(1);
            else
                %display('LB=3, LD=0: Not fixed');
                REM1=[]; REM2=[];
                break;
            end
        elseif LB>3 && LD==0            
            %display('LB>3 && LD==0');
            REM1=[]; REM2=[];
            break;
        end        

        x1=zeros(1,LD);  y1=zeros(1,LD);  a1=zeros(1,LD);  t1=zeros(1,LD);  p1=zeros(LD,6);    
        x2=zeros(1,LD);  y2=zeros(1,LD);  a2=zeros(1,LD);  t2=zeros(1,LD);  p2=zeros(LD,6);    
        x3=zeros(1,LD);  y3=zeros(1,LD);  a3=zeros(1,LD);  t3=zeros(1,LD);  p3=zeros(LD,6);        

        for i=1:LD

             x1(i)=XO(NumInter(i)); y1(i)=YO(NumInter(i)); 
             a1(i)=AO(NumInter(i)); t1(i)=TO(NumInter(i)); 
             p1(i,:)=PO(NumInter(i),:);

             if NumInter(i)==LB
                x2(i)=XO(1); y2(i)=YO(1); 
                a2(i)=AO(1); t2(i)=TO(1);  
                p2(i,:)=PO(1,:);
             else
                x2(i)=XO(NumInter(i)+1); y2(i)=YO(NumInter(i)+1); 
                a2(i)=AO(NumInter(i)+1); t2(i)=TO(NumInter(i)+1);            
                p2(i,:)=PO(NumInter(i)+1,:);
            end

             x3(i)=Xinter(i);
             y3(i)=Yinter(i);
             a3(i)=Ainter(i);
             t3(i)=Tinter(i);
             p3(i,:)=Pinter(i,:);

             
             if t1(i)==1
                 rem=rem+1;
                 rem1=rem1+1;
                 %line([x1(i) x3(i)],[y1(i) y3(i)],'Color',[.2 .2 .2],'LineWidth',1,'LineStyle','-');
                 REM1(rem1,:)=[x1(i) y1(i) x3(i) y3(i)];
             elseif t1(i)==2
                 rem=rem+1;
                 rem1=rem1+1;
                 %line([x1(i) x3(i)],[y1(i) y3(i)],'Color',[.2 .2 .2],'LineWidth',1,'LineStyle','-');
                 REM1(rem1,:)=[x1(i) y1(i) x3(i) y3(i)];
             elseif t1(i)==-1
                 rem=rem+1;
                 rem2=rem2+1;
                 %line([x1(i) x3(i)],[y1(i) y3(i)],'Color',[.8 .8 .8],'LineWidth',1,'LineStyle','-');         
                 REM2(rem2,:)=[x1(i) y1(i) x3(i) y3(i)];
             elseif t1(i)==0
                 rem=rem+1;
                 rem1=rem1+1;
                 %ParabPlotX(x1(i),y1(i),p1(i,:),x3(i),y3(i),XB,YB,AL,AR,[.2 .2 .2]);              
                 REM1(rem1,:)=[x1(i) y1(i) x3(i) y3(i)];             
             end

             if t2(i)==1
                rem=rem+1;
                rem1=rem1+1;
                %line([x2(i) x3(i)],[y2(i) y3(i)],'Color',[.2 .2 .2],'LineWidth',1,'LineStyle','-');
                REM1(rem1,:)=[x2(i) y2(i) x3(i) y3(i)];
             elseif t2(i)==2         
                rem=rem+1;
                rem1=rem1+1;
                %line([x2(i) x3(i)],[y2(i) y3(i)],'Color',[.2 .2 .2],'LineWidth',1,'LineStyle','-');
                REM1(rem1,:)=[x2(i) y2(i) x3(i) y3(i)];
             elseif t2(i)==-1  
                rem=rem+1;
                rem2=rem2+1;
                %line([x2(i) x3(i)],[y2(i) y3(i)],'Color',[.8 .8 .8],'LineWidth',1,'LineStyle','-');         
                REM2(rem2,:)=[x2(i) y2(i) x3(i) y3(i)];
             elseif t2(i)==0
                rem=rem+1;
                rem1=rem1+1;
                %ParabPlotX(x2(i),y2(i),p2(i,:),x3(i),y3(i),XB,YB,AL,AR,[.2 .2 .2]);
                REM1(rem1,:)=[x2(i) y2(i) x3(i) y3(i)];
             end         
        end        

        XOn=[]; YOn=[]; AOn=[]; TOn=[]; POn=[]; NOn=[]; 
        j=0; q=0;
        for i=1:LB
            n1=i;
            if i==1
                n2=LB;
            else
                n2=i-1;
            end
            indx1=sum(n1==NumInter);
            indx2=sum(n2==NumInter);
            if indx1==0 && indx2==0
                q=q+1;
                XOn=[XOn XO(n1)];   YOn=[YOn YO(n1)];   
                AOn=[AOn AO(n1)];   TOn=[TOn TO(n1)];
                POn=[POn; PO(n1,:)];            

            elseif indx1==1 && indx2==0
                j=j+1;
                q=q+1;
                XOn=[XOn x3(j)];    YOn=[YOn y3(j)];   
                AOn=[AOn a3(j)];    TOn=[TOn t3(j)];
                POn=[POn; p3(j,:)];            
                NOn=[NOn q];                
            end        
        end
        
        if ~isempty(NOn)
            MEM=NOn;
            MEMm=NOn-1;
            MEMp=NOn+1;
            MEMm(MEMm==0)=length(XOn);
            MEMp(MEMp==(length(XOn)+1))=1;
            MEMm=[MEMm(2:end) MEMm(1)];
            MEMp=[MEMp(end) MEMp(1:(end-1))];
            MEM=sort([MEM MEMm(MEMm~=MEM) MEMp(MEMp~=MEM)]);
        end
        
        XO=XOn; YO=YOn; AO=AOn; TO=TOn; PO=POn;
        LB=length(XO); 
        NO=NOn;
        
    end

    if LB==2
        %- need to add chacking if rem1 or rem2 !!!!!!!!!!!!!!!!!!!!!!!
        rem=rem+1;
        rem1=rem1+1;
        %line(XO,YO,'Color',[.2 .2 .2],'LineWidth',1,'LineStyle','-'); 
        REM1(rem1,:)=[XO(1) YO(1) XO(2) YO(2)];
    end
    
    Nedge=2*(Nconvex+Nconcave)+2*Nconcave-3;
    Nvert=2*(Nconvex+Nconcave)+1*Nconcave-2;
    
    if LB<3 && rem==Nedge
        REM=[REM1; REM2];
        RR=[REM(:,1:2); REM(:,3:4)];
        R=RR(:,1).^2+RR(:,2).^2;
        U=unique(R);
        if size(U,1)==Nvert
            Status=1;            
        else
            REM1=[];  REM2=[];
            %disp('LB & rem & Nedge');
            %disp([LB rem Nedge]);
            %disp('size(U,1) & Nvert');
            %disp([size(U,1) Nvert]);
        end
    else
        %disp('LB & rem & Nedge');
        %disp([LB rem Nedge]);
        REM1=[];  REM2=[];
    end
    
    
    cla reset;    
    set(handles.progress1,'Box','on','XTick',[],'YTick',[],'Visible','on');
    %cla(handles.progress1);
    %set(handles.progress1,'Visible','off');
    
    
