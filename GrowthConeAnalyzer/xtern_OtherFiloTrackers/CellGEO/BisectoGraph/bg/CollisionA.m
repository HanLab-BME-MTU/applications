function [NumInter,Xinter,Yinter,Ainter,Tinter,Pinter]=CollisionA(XB,YB,MB,AL,AR,XO,YO,AO,TO,PO,MEM)

NumInter=[];    Xinter=[];    Yinter=[];    Ainter=[];    Tinter=[];    Pinter=[];    %dMEM=[];   

LB=length(XO);

dc=[]; ic=[];
xc=[]; yc=[]; ac=[]; tc=[]; pc=[];
LoR=[];  % 0 - left (first), 1 - right (second)

LN=length(MEM);
for ii=1:LN
    
    i=MEM(ii);
    if i==1
        xo1=XO(LB); xo2=XO(i); xo3=XO(i+1);
        yo1=YO(LB); yo2=YO(i); yo3=YO(i+1);        
        ao1=AO(LB); ao2=AO(i); ao3=AO(i+1);        
        to1=TO(LB); to2=TO(i); to3=TO(i+1);
        xm1=PO(LB,1); xm2=PO(i,1); xm3=PO(i+1,1);
        ym1=PO(LB,2); ym2=PO(i,2); ym3=PO(i+1,2);
        prm1=PO(LB,:); prm2=PO(i,:); prm3=PO(i+1,:);
        ko1=tan(AO(LB)); 
        ko2=tan(AO(i)); 
        ko3=tan(AO(i+1));                
    elseif i==LB
        xo1=XO(i-1); xo2=XO(i); xo3=XO(1);
        yo1=YO(i-1); yo2=YO(i); yo3=YO(1);
        ao1=AO(i-1); ao2=AO(i); ao3=AO(1);        
        to1=TO(i-1); to2=TO(i); to3=TO(1);
        xm1=PO(i-1,1); xm2=PO(i,1); xm3=PO(1,1);
        ym1=PO(i-1,2); ym2=PO(i,2); ym3=PO(1,2);
        prm1=PO(i-1,:); prm2=PO(i,:); prm3=PO(1,:);
        ko1=tan(AO(i-1)); 
        ko2=tan(AO(i)); 
        ko3=tan(AO(1));
    else
        xo1=XO(i-1); xo2=XO(i); xo3=XO(i+1);
        yo1=YO(i-1); yo2=YO(i); yo3=YO(i+1);
        ao1=AO(i-1); ao2=AO(i); ao3=AO(i+1);        
        to1=TO(i-1); to2=TO(i); to3=TO(i+1);
        xm1=PO(i-1,1); xm2=PO(i,1); xm3=PO(i+1,1);
        ym1=PO(i-1,2); ym2=PO(i,2); ym3=PO(i+1,2);
        prm1=PO(i-1,:); prm2=PO(i,:); prm3=PO(i+1,:);
        ko1=tan(AO(i-1)); 
        ko2=tan(AO(i)); 
        ko3=tan(AO(i+1));
    end 
    
    %--------------------------------- f ----------------------------------
    
    if ii==1 || (MEM(ii)-MEM(ii-1)>1)
        if     to1==-1 && to2==-1
            indxf=0;
        elseif to1== 0 && to2== 0        
            [df,dt,xf,yf,af,tf,pf,status] = ParabParabLine(xo2,yo2,ao2,prm2,xo1,yo1,ao1,prm1,XB,YB,AL,AR);
            if status==1
                indxf=inpolygon(xf,yf,[XB'; XB(1)],[YB'; YB(1)]);
            else
                indxf=0;
            end          
        elseif to1== 0 && to2== 1
            [df,xf,yf,af,tf,pf,status] = LineParabParabA(xm2,ym2,ao2,prm2,xo1,yo1,ao1,prm1,XB,YB,AL,AR);
            df=sqrt((xf-xo2)^2+(yf-yo2)^2);
            if status==1
                indxf=inpolygon(xf,yf,[XB'; XB(1)],[YB'; YB(1)]);
            else
                indxf=0;
            end     
        elseif to1== 1 && to2== 0
            [df,xf,yf,af,tf,pf,status] = LineParabParabA(xm1,ym1,ao1,prm1,xo2,yo2,ao2,prm2,XB,YB,AL,AR);        
            if status==1
                indxf=inpolygon(xf,yf,[XB'; XB(1)],[YB'; YB(1)]);
            else
                indxf=0;
            end
        elseif to1== 0 && to2== 2
            [df,xf,yf,af,tf,pf,status] = MLineParabParabA(xm2,ym2,ao2,prm2,xo1,yo1,ao1,prm1,XB,YB,AL,AR);
            df=sqrt((xf-xo2)^2+(yf-yo2)^2);
            if status==1
                indxf=inpolygon(xf,yf,[XB'; XB(1)],[YB'; YB(1)]);
            else
                indxf=0;
            end                  
        elseif to1== 2 && to2== 0
            [df,xf,yf,af,tf,pf,status] = MLineParabParabA(xm1,ym1,ao1,prm1,xo2,yo2,ao2,prm2,XB,YB,AL,AR);        
            if status==1
                indxf=inpolygon(xf,yf,[XB'; XB(1)],[YB'; YB(1)]);
            else
                indxf=0;
            end               
        elseif to1== 0 && to2==-1
            [df,xf,yf,af,tf,pf,status] = LineParabLine(xo1,yo1,ao1,to1,prm1,xo2,yo2,ao2,to2,prm2,XB,YB,AL,AR);
            df=sqrt((xf-xo2)^2+(yf-yo2)^2);
            if status==1
                indxf=inpolygon(xf,yf,[XB'; XB(1)],[YB'; YB(1)]);
            else
                indxf=0;
            end               
        elseif  to1==-1 && to2==0
            [df,xf,yf,af,tf,pf,status] = LineParabLine(xo1,yo1,ao1,to1,prm1,xo2,yo2,ao2,to2,prm2,XB,YB,AL,AR);
            if status==1
                indxf=inpolygon(xf,yf,[XB'; XB(1)],[YB'; YB(1)]);
            else
                indxf=0;
            end        
        else 
            xf=(ko2*xo2-ko1*xo1)/(ko2-ko1)-(yo2-yo1)/(ko2-ko1);
            yf=ko2*ko1*(xo2-xo1)/(ko2-ko1)-(ko1*yo2-ko2*yo1)/(ko2-ko1);

            if     to1== 1 && to2==-1
                af=ao1;
                tf=0;
                pf=LineLineParab(xo1,yo1,ao1,to1,prm1,xo2,yo2,ao2,to2,prm2,XB,YB,AL,AR);
            elseif to1== 2 && to2==-1
                af=ao1;
                tf=0;
                pf=MLineLineParab(xo1,yo1,ao1,to1,prm1,xo2,yo2,ao2,to2,prm2,XB,YB,AL,AR);
                if isempty(pf)
                    NumInter=-1;
                    return;
                end
            elseif to1==-1 && to2== 1
                af=ao2;
                tf=0;
                pf=LineLineParab(xo1,yo1,ao1,to1,prm1,xo2,yo2,ao2,to2,prm2,XB,YB,AL,AR);            
            elseif to1==-1 && to2== 2
                af=ao2;
                tf=0;
                pf=MLineLineParab(xo1,yo1,ao1,to1,prm1,xo2,yo2,ao2,to2,prm2,XB,YB,AL,AR);
                if isempty(pf)
                    NumInter=-1;
                    return;
                end
            elseif to1== 1 && to2== 1            
                tf=1;
                [af,pf] = LineLineLine(ao1,prm1,ao2,prm2,XB,YB,AL,AR);
            elseif to1== 1 && to2== 2            
                tf=1;
                [af,pf] = MLineLineLine(ao1,prm1,ao2,prm2); 
            elseif to1== 2 && to2== 1            
                tf=1;
                [af,pf] = MLineLineLine(ao1,prm1,ao2,prm2);
            elseif to1== 2 && to2== 2            
                tf=2;
                [af,pf] = MMLineLineLine(ao1,prm1,ao2,prm2,XB,YB,AL,AR);
                if isempty(pf)
                    NumInter=-1;
                    return;
                end
            end

            df=sqrt((xf-xo2)^2+(yf-yo2)^2);
            an=atan2(yf-yo2,xf-xo2);
            if sign(cos(an-ao2))==-1
                indxf=0;
            else
                indxf=inpolygon(xf,yf,[XB'; XB(1)],[YB'; YB(1)]);
            end

        end

        if indxf==1

            check = LinearityCheck(xo2,yo2,to2,prm2,xf,yf,XB,YB,MB);
            if check==0
                indxf=0;
            end
        end 
    else
        
        indxf=indxs;
        if indxf
            df=dm; xf=xs; yf=ys;
            af=as; tf=ts; pf=ps;        
        end
    end
    %--------------------------------- s ----------------------------------
    
    if     to2==-1 && to3==-1  
        indxs=0;
    elseif to2== 0 && to3== 0 
        [ds,dm,xs,ys,as,ts,ps,status] = ParabParabLine(xo2,yo2,ao2,prm2,xo3,yo3,ao3,prm3,XB,YB,AL,AR);        
        if status==1
            indxs=inpolygon(xs,ys,[XB'; XB(1)],[YB'; YB(1)]);
        else
            indxs=0;
        end      
    elseif to2== 0 && to3== 1
        [ds,xs,ys,as,ts,ps,status] = LineParabParabA(xm3,ym3,ao3,prm3,xo2,yo2,ao2,prm2,XB,YB,AL,AR);        
        dm=sqrt((xs-xo3)^2+(ys-yo3)^2);
        if status==1
            indxs=inpolygon(xs,ys,[XB'; XB(1)],[YB'; YB(1)]);
        else
            indxs=0;
        end     
    elseif to2== 1 && to3== 0
        [dm,xs,ys,as,ts,ps,status] = LineParabParabA(xm2,ym2,ao2,prm2,xo3,yo3,ao3,prm3,XB,YB,AL,AR);
        ds=sqrt((xs-xo2)^2+(ys-yo2)^2);
        if status==1
            indxs=inpolygon(xs,ys,[XB'; XB(1)],[YB'; YB(1)]);
        else
            indxs=0;
        end
    elseif to2== 0 && to3== 2
        [ds,xs,ys,as,ts,ps,status] = MLineParabParabA(xm3,ym3,ao3,prm3,xo2,yo2,ao2,prm2,XB,YB,AL,AR);        
        dm=sqrt((xs-xo3)^2+(ys-yo3)^2);
        if status==1
            indxs=inpolygon(xs,ys,[XB'; XB(1)],[YB'; YB(1)]);
        else
            indxs=0;
        end         
    elseif to2== 2 && to3== 0
        [dm,xs,ys,as,ts,ps,status] = MLineParabParabA(xm2,ym2,ao2,prm2,xo3,yo3,ao3,prm3,XB,YB,AL,AR);
        ds=sqrt((xs-xo2)^2+(ys-yo2)^2);
        if status==1
            indxs=inpolygon(xs,ys,[XB'; XB(1)],[YB'; YB(1)]);
        else
            indxs=0;
        end        
    elseif to2== 0 && to3==-1
        [ds,xs,ys,as,ts,ps,status] = LineParabLine(xo2,yo2,ao2,to2,prm2,xo3,yo3,ao3,to3,prm3,XB,YB,AL,AR);
        dm=sqrt((xs-xo3)^2+(ys-yo3)^2);
        if status==1
            indxs=inpolygon(xs,ys,[XB'; XB(1)],[YB'; YB(1)]);
        else
            indxs=0;
        end
    elseif to2==-1 && to3== 0
        [dm,xs,ys,as,ts,ps,status] = LineParabLine(xo2,yo2,ao2,to2,prm2,xo3,yo3,ao3,to3,prm3,XB,YB,AL,AR);
        ds=sqrt((xs-xo2)^2+(ys-yo2)^2);
        if status==1
            indxs=inpolygon(xs,ys,[XB'; XB(1)],[YB'; YB(1)]);
        else
            indxs=0;
        end
    else 
        xs=(ko2*xo2-ko3*xo3)/(ko2-ko3)-(yo2-yo3)/(ko2-ko3);
        ys=ko2*ko3*(xo2-xo3)/(ko2-ko3)-(ko3*yo2-ko2*yo3)/(ko2-ko3);
        
        if     to2== 1 && to3==-1
            as=ao2;
            ts=0;
            ps=LineLineParab(xo2,yo2,ao2,to2,prm2,xo3,yo3,ao3,to3,prm3,XB,YB,AL,AR);
        elseif to2== 2 && to3==-1
            as=ao2;
            ts=0;
            ps=MLineLineParab(xo2,yo2,ao2,to2,prm2,xo3,yo3,ao3,to3,prm3,XB,YB,AL,AR);  
            if isempty(ps)
                NumInter=-1;
                return;
            end
        elseif to2==-1 && to3== 1
            as=ao3;
            ts=0;
            ps=LineLineParab(xo2,yo2,ao2,to2,prm2,xo3,yo3,ao3,to3,prm3,XB,YB,AL,AR);
        elseif to2==-1 && to3== 2
            as=ao3;
            ts=0;
            ps=MLineLineParab(xo2,yo2,ao2,to2,prm2,xo3,yo3,ao3,to3,prm3,XB,YB,AL,AR);   
            if isempty(ps)
                NumInter=-1;
                return;
            end
        elseif to2== 1 && to3== 1                   
            ts=1;
            [as,ps] = LineLineLine(ao2,prm2,ao3,prm3,XB,YB,AL,AR);
        elseif to2== 1 && to3== 2                   
            ts=1;
            [as,ps] = MLineLineLine(ao3,prm3,ao2,prm2);
        elseif to2== 2 && to3== 1                   
            ts=1;
            [as,ps] = MLineLineLine(ao3,prm3,ao2,prm2);
        elseif to2== 2 && to3== 2                   
            ts=2;
            [as,ps] = MMLineLineLine(ao3,prm3,ao2,prm2,XB,YB,AL,AR);
            if isempty(ps)
                NumInter=-1;
                return;
            end
        end
        
        ds=sqrt((xs-xo2)^2+(ys-yo2)^2);
        dm=sqrt((xs-xo3)^2+(ys-yo3)^2);
        an=atan2(ys-yo2,xs-xo2);
        if sign(cos(an-ao2))==-1
            indxs=0;
        else
            indxs=inpolygon(xs,ys,[XB'; XB(1)],[YB'; YB(1)]);
        end
        
    end
    
    if indxs==1
        check = LinearityCheck(xo2,yo2,to2,prm2,xs,ys,XB,YB,MB);
        if check==0
            indxs=0;
        end        
    end  
    
    %----------------------------------------------------------------------
    
    if indxf && indxs
        ic=[ic i];
        if df<ds
           dc=[dc df];  xc=[xc xf];  yc=[yc yf];  
           ac=[ac af];  tc=[tc tf];  pc=[pc; pf]; LoR=[LoR 0];
        else
           dc=[dc ds];  xc=[xc xs];  yc=[yc ys];
           ac=[ac as];  tc=[tc ts];  pc=[pc; ps]; LoR=[LoR 1];
        end              
    elseif indxf && ~indxs
        ic=[ic i];  dc=[dc df];  xc=[xc xf];  yc=[yc yf];  
                    ac=[ac af];  tc=[tc tf];  pc=[pc; pf]; LoR=[LoR 0];
    elseif ~indxf && indxs
        ic=[ic i];  dc=[dc ds];  xc=[xc xs];  yc=[yc ys];  
                    ac=[ac as];  tc=[tc ts];  pc=[pc; ps]; LoR=[LoR 1];
    end
    
    %----------------------------------------------------------------------    
    
end

LD=length(ic);

if LD==0
    NumInter=[];        
elseif LD==1
    NumInter=[NumInter 0];   
%{    
elseif LD==2
    if ic(1)==1 && ic(2)~=2 %!!!!!!!!!!!!!!!!!!!!!        
        n1=2; n2=1;
    else
        n1=1; n2=2;
    end
    dis=sqrt((xc(n1)-xc(n2))^2+(yc(n1)-yc(n2))^2);
    if dis<1e-9 %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        NumInter=[NumInter ic(n1)];
        Xinter=[Xinter xc(n1)];
        Yinter=[Yinter yc(n1)]; 
        Ainter=[Ainter ac(n1)]; 
        Tinter=[Tinter tc(n1)]; 
        Pinter=[Pinter; pc(n1,:)]; 
    end
%}
else
    for j=1:LD
        
        edd=0;
        if j==LD
           if ic(LD)==LB && ic(1)==1 && LoR(LD)==1 && LoR(1)==0
               edd=1;
           end
        else
            if ic(j+1)-ic(j)==1 && LoR(j)==1 && LoR(j+1)==0
                edd=1;
            end
        end
        
        if edd
            NumInter=[NumInter ic(j)];
            Xinter=[Xinter xc(j)];
            Yinter=[Yinter yc(j)]; 
            Ainter=[Ainter ac(j)]; 
            Tinter=[Tinter tc(j)]; 
            Pinter=[Pinter; pc(j,:)]; 
        end        
    end
end     