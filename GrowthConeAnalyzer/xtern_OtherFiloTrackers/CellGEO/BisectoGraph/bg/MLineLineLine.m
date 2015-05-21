function [ac,pc] = MLineLineLine(a1,prm1,a2,prm2)
 
    ac=[]; pc=[];
    
    xq1=prm1(3); yq1=prm1(4);
    xq2=prm2(3); yq2=prm2(4);
    
    if abs(xq1)+abs(yq1)>1e-9 && abs(xq2)+abs(yq2)>1e-9
        M=[[-sin(a1-pi/2) , cos(a1-pi/2)      ];...
           [-sin(a2-pi/2) , cos(a2-pi/2) ]];
        m=[ -sin(a1-pi/2)*xq1 + cos(a1-pi/2)*yq1  ;...
            -sin(a2-pi/2)*xq2 + cos(a2-pi/2)*yq2];

        DM=det(M);
        if abs(DM)<1e-12
             return;
        end
        temp=M\m;
        xm=temp(1);
        ym=temp(2);
        
        xp1=2*xq1-xm; yp1=2*yq1-ym;
        xp2=2*xq2-xm; yp2=2*yq2-ym;
        xq3=(xp1+xp2)/2; yq3=(yp1+yp2)/2;
        pc=[xq3 yq3 xq3 yq3 0 0];
        ac=-pi/2+atan2(yp2-yp1,xp2-xp1);
        if ac<0
            ac=ac+2*pi;
        end
        
        if cos(ac-a1)<0 || cos(ac-a2)<0
            ac=ac-pi;            
        end
        if ac<0
            ac=ac+2*pi;            
        end  
        return;        
    end

    xm1=prm1(1); ym1=prm1(2);
    xm2=prm2(1); ym2=prm2(2);
    
    gam=atan2(ym2-ym1,xm2-xm1);
    if gam<0
        gam=gam+2*pi;
    end
    ac=a1+a2-gam-pi/2;
    if ac<0
        ac=ac+2*pi;
    elseif ac>(2*pi)
        ac=ac-2*pi;
    end
    phi1=2*a1-gam;
    if phi1<0
        phi1=phi1+2*pi;
    elseif phi1>(2*pi)
        phi1=phi1-2*pi;
    end    
    phi2=2*a2-gam-pi;
    if phi2<0
        phi2=phi2+2*pi;
    elseif phi2>(2*pi)
        phi2=phi2-2*pi;
    end

    am=atan2(0.5*sin(a1)+0.5*sin(a2),0.5*cos(a1)+0.5*cos(a2));
    if cos(ac-am)<0
        ac=ac-pi;
    end
    if ac<0
        ac=ac+2*pi;            
    end    

    M=[[-sin(phi1) , cos(phi1)      ];...
       [-sin(phi2) , cos(phi2) ]];
    m=[ -sin(phi1)*xm1 + cos(phi1)*ym1  ;...
        -sin(phi2)*xm2 + cos(phi2)*ym2];

    DM=det(M);
    if abs(DM)<1e-12
        return;
    end
    temp=M\m;
    xt=temp(1);
    yt=temp(2);
    pc=[xt yt 0 0 0 0];


return;

