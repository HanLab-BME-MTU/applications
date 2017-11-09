function [REM,Status]=FindGraph(Xb,Yb,Lb)

    XB=Xb; YB=Yb; LB=Lb;
    MB=[max(XB) min(XB) max(YB) min(YB)];

    AL=zeros(1,LB); AR=zeros(1,LB);
    XO=[]; YO=[]; AO=[]; TO=[]; PO=[];

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
            return;
        elseif pntyp==1;
            gam=atan2(sin(AL(i))+sin(AR(i)),cos(AL(i))+cos(AR(i)));
            if gam<0
                gam=gam+2*pi;
            end
            XO=[XO x1];  YO=[YO y1];
            AO=[AO gam]; TO=[TO pntyp];        
            PO=[PO; [x1 y1 0 0 i i]]; 
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
        end   

    end
    LB=length(XO);
    NO=1:LB;
    MEM=NO;
        
    %--------------------------------------------------------------------------
    REM=zeros(2*LB-3,4);
    rem=0;
    e=0;    
    while LB>2    
        e=e+1;
        [NumInter,Xinter,Yinter,Ainter,Tinter,Pinter]=CollisionA(XB,YB,MB,AL,AR,XO,YO,AO,TO,PO,MEM);
       
        if isempty(NumInter) || sum(NumInter==0)>0 
            MEM=1:length(XO);
            [NumInter,Xinter,Yinter,Ainter,Tinter,Pinter]=CollisionA(XB,YB,MB,AL,AR,XO,YO,AO,TO,PO,MEM);
            if isempty(NumInter) 
                display('NumInter is empty');
                break;
            elseif NumInter==0
                display('NumInter is zero');
                break;
            end
        elseif NumInter==-1
            break;
        end
        
        LD=length(NumInter);

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
                 line([x1(i) x3(i)],[y1(i) y3(i)],'Color',[.2 .2 .2],'LineWidth',1,'LineStyle','-');
                 REM(rem,:)=[x1(i) y1(i) x3(i) y3(i)];
             elseif t1(i)==2
                 rem=rem+1;
                 line([x1(i) x3(i)],[y1(i) y3(i)],'Color',[.2 .2 .2],'LineWidth',1,'LineStyle','-');
                 REM(rem,:)=[x1(i) y1(i) x3(i) y3(i)];
             elseif t1(i)==-1
                 line([x1(i) x3(i)],[y1(i) y3(i)],'Color',[.8 .8 .8],'LineWidth',1,'LineStyle','-');                          
             elseif t1(i)==0
                 rem=rem+1;
                 ParabPlotX(x1(i),y1(i),p1(i,:),x3(i),y3(i),XB,YB,AL,AR,[.2 .2 .2]);              
                 REM(rem,:)=[x1(i) y1(i) x3(i) y3(i)];             
             end

             if t2(i)==1
                rem=rem+1;
                line([x2(i) x3(i)],[y2(i) y3(i)],'Color',[.2 .2 .2],'LineWidth',1,'LineStyle','-');
                REM(rem,:)=[x2(i) y2(i) x3(i) y3(i)];
             elseif t2(i)==2         
                rem=rem+1; 
                line([x2(i) x3(i)],[y2(i) y3(i)],'Color',[.2 .2 .2],'LineWidth',1,'LineStyle','-');
                REM(rem,:)=[x2(i) y2(i) x3(i) y3(i)];
             elseif t2(i)==-1         
                line([x2(i) x3(i)],[y2(i) y3(i)],'Color',[.8 .8 .8],'LineWidth',1,'LineStyle','-');                         
             elseif t2(i)==0
                rem=rem+1;
                ParabPlotX(x2(i),y2(i),p2(i,:),x3(i),y3(i),XB,YB,AL,AR,[.2 .2 .2]);
                REM(rem,:)=[x2(i) y2(i) x3(i) y3(i)];
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
            MEM=unique([MEM MEMm MEMp]);
        end        
       
        XO=XOn; YO=YOn; AO=AOn; TO=TOn; PO=POn;
        LB=length(XO); 
        NO=NOn;        
    end

    if LB==2
        rem=rem+1;
        line(XO,YO,'Color',[.2 .2 .2],'LineWidth',1,'LineStyle','-'); 
        REM(rem,:)=[XO(1) YO(1) XO(2) YO(2)];
    end
    
    if LB<3
        Status=1;
    else
        Status=0;
    end

    im=find(sum(REM,2)==0,1)-1;
    if isempty(im)
        im=size(REM,1);
    end
    REM=REM(1:im,:);
    
           
        
    