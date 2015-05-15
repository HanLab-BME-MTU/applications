function [dc,dm,xc,yc,ac,tc,pc,status] = ParabParabLine(x1,y1,a1,prm1,x2,y2,a2,prm2,XB,YB,AL,AR)

    pi=3.1415926535897932384626433832795;   

    dc=[]; dm=[]; xc=[]; yc=[]; ac=[]; tc=[]; pc=[]; status=0; status1=0; status2=0;
    
    ksi1=prm1(3); K1=prm1(5); N1=prm1(6);
    ksi2=prm2(3); K2=prm2(5); N2=prm2(6);
    
    XD1=XB(N1); YD1=YB(N1);
    XD2=XB(N2); YD2=YB(N2);
    XF1=XB(K1); YF1=YB(K1);
    XF2=XB(K2); YF2=YB(K2);
    
    if abs(cos(ksi1-AL(N1)))<1e-9 && abs(cos(ksi1-AR(N1)))>1e-9
        AD1=AL(N1);
    elseif abs(cos(ksi1-AL(N1)))>1e-9 && abs(cos(ksi1-AR(N1)))<1e-9
        AD1=AR(N1);
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else
        return;
    end

    if abs(cos(ksi2-AL(N2)))<1e-9 && abs(cos(ksi2-AR(N2)))>1e-9
        AD2=AL(N2);
    elseif abs(cos(ksi2-AL(N2)))>1e-9 && abs(cos(ksi2-AR(N2)))<1e-9
        AD2=AR(N2);
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else
        return;
    end
    
    M=[[-sin(AD1), cos(AD1)];...
       [-sin(AD1-pi/2), cos(AD1-pi/2)]];
    m=[ -sin(AD1)*XD1 + cos(AD1)*YD1;...
        -sin(AD1-pi/2)*XF1 + cos(AD1-pi/2)*YF1];

    DM=det(M);
    if abs(DM)<1e-12
        return;
    end
    temp=M\m;
    xs1=temp(1);
    ys1=temp(2);

    M=[[-sin(AD2), cos(AD2)];...
       [-sin(AD2-pi/2), cos(AD2-pi/2)]];
    m=[ -sin(AD2)*XD2 + cos(AD2)*YD2;...
        -sin(AD2-pi/2)*XF2 + cos(AD2-pi/2)*YF2];

    DM=det(M);
    if abs(DM)<1e-12
        return;
    end
    temp=M\m;
    xs2=temp(1);
    ys2=temp(2);

    %XO1=(XF1+xs1)/2; YO1=(YF1+ys1)/2;
    %XO2=(XF2+xs2)/2; YO2=(YF2+ys2)/2;
    P1=sqrt((XF1-xs1)^2+(YF1-ys1)^2);
    P2=sqrt((XF2-xs2)^2+(YF2-ys2)^2);   
    
    M=[[-sin(AD1-pi/2) , cos(AD1-pi/2)  ];...
       [-sin(AD1) , cos(AD1) ]];
    m=[ -sin(AD1-pi/2)*XF1 + cos(AD1-pi/2)*YF1  ;...
        -sin(AD1)*x1 + cos(AD1)*y1 ];

    DM=det(M);
    if abs(DM)<1e-12
        return;
    end
    temp=M\m;
    xt1=temp(1);
    yt1=temp(2);

    sp1=sign((x1-xt1)*sin(ksi1) - (y1-yt1)*cos(ksi1));
    dp1=sp1*sqrt((x1-xt1)^2+(y1-yt1)^2); 
    
    M=[[-sin(AD2-pi/2) , cos(AD2-pi/2)  ];...
       [-sin(AD2) , cos(AD2) ]];
    m=[ -sin(AD2-pi/2)*XF2 + cos(AD2-pi/2)*YF2  ;...
        -sin(AD2)*x2 + cos(AD2)*y2 ];

    DM=det(M);
    if abs(DM)<1e-12
        return;
    end
    temp=M\m;
    xt2=temp(1);
    yt2=temp(2);

    sp2=sign((x2-xt2)*sin(ksi2) - (y2-yt2)*cos(ksi2));
    dp2=sp2*sqrt((x2-xt2)^2+(y2-yt2)^2); 

    
    if K1==K2
        
        ac=atan2(sin(AD1)+sin(AD2),cos(AD1)+cos(AD2));
        if ac<0
            ac=ac+2*pi;
        end

        if sign(sin(ac-AD1))==1 && sign(sin(ac-AD2))==-1
            N3=N2; M3=N1;
        elseif sign(sin(ac-AD1))==-1 && sign(sin(ac-AD2))==1
            N3=N1; M3=N2;
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else
            return;
        end 
        
        M=[[-sin(ac-pi/2) , cos(ac-pi/2)  ];...
           [-sin(AD2) , cos(AD2) ]];
        m=[ -sin(ac-pi/2)*XD1 + cos(ac-pi/2)*YD1  ;...
            -sin(AD2)*XD2 + cos(AD2)*YD2 ];

        DM=det(M);
        if abs(DM)<1e-12
            return;
        end
        temp=M\m;
        xc1=(XD1+temp(1))/2;
        yc1=(YD1+temp(2))/2;
        
        
        [XC,YC,DC,stfd]=LineParabInterA(XF1,YF1,AD1,XD1,YD1,ac,xc1,yc1);
        
        xi=[]; yi=[]; dd=[]; dw=[];
        for i=1:2
            if stfd
                
                xc=XC(i); yc=YC(i); dc1=DC(i);
                
                M=[[-sin(AD2-pi/2) , cos(AD2-pi/2)  ];...
                   [-sin(AD2) , cos(AD2) ]];
                m=[ -sin(AD2-pi/2)*XF2 + cos(AD2-pi/2)*YF2  ;...
                    -sin(AD2)*xc + cos(AD2)*yc ];

                DM=det(M);
                if abs(DM)<1e-12
                    return;
                end
                temp=M\m;
                xt=temp(1);
                yt=temp(2);

                sc2=sign((xc-xt)*sin(ksi2) - (yc-yt)*cos(ksi2));
                dc2=sc2*sqrt((xc-xt)^2+(yc-yt)^2);                 
                
                
                if cos(ksi1-a1)<0 && dp1>0
                    if dc1<dp1
                        status1=1;
                    else
                        status1=0;
                    end
                elseif cos(ksi1-a1)>0 && dp1>0
                    if dc1>dp1
                        status1=1;
                    else
                        status1=0;
                    end
                elseif cos(ksi1-a1)<0 && dp1<0    
                    if dc1>dp1
                        status1=1;
                    else
                        status1=0;
                    end
                elseif cos(ksi1-a1)>0 && dp1<0  
                    if dc1<dp1
                        status1=1;
                    else
                        status1=0;
                    end
                end
                
                if cos(ksi2-a2)<0 && dp2>0
                    if dc2<dp2
                        status2=1;
                    else
                        status2=0;
                    end
                elseif cos(ksi2-a2)>0 && dp2>0
                    if dc2>dp2
                        status2=1;
                    else
                        status2=0;
                    end
                elseif cos(ksi2-a2)<0 && dp2<0    
                    if dc2>dp2
                        status2=1;
                    else
                        status2=0;
                    end
                elseif cos(ksi2-a2)>0 && dp2<0  
                    if dc2<dp2
                        status2=1;
                    else
                        status2=0;
                    end
                end
                
                if status1==1 && status2==1
                    xi=[xi xc];
                    yi=[yi yc];
                    dd=[dd abs(dc1-dp1)];
                    dw=[dw abs(dc2-dp2)];
                end
            end
            
        end
        
        if isempty(dd)               
            return;
        else
            M=[[-sin(AD1), cos(AD1)];...
               [-sin(AD2), cos(AD2)]];
            m=[ -sin(AD1)*XD1 + cos(AD1)*YD1;...
                -sin(AD2)*XD2 + cos(AD2)*YD2];

            DM=det(M);
            if abs(DM)<1e-12
                return;
            end
            temp=M\m;
            xm=temp(1);
            ym=temp(2);
            [dc,it]=min(dd);
            dm=dw(it);
            xc=xi(it);
            yc=yi(it);            
            tc=1;
            pc=[xm ym 0 0 N3 M3];                                       
            
            status=1;            
        end
        
        
        
    else
        
        if abs(sin(AD1-AD2))>1e-9      
            return;
        end

        xq=(XF1+XF2)/2;
        yq=(YF1+YF2)/2;
        ac=-pi/2+atan2(YF2-YF1,XF2-XF1);
        if ac<0
            ac=ac+2*pi;
        end
        
        d=sqrt((xs1-xs2)^2+(ys1-ys2)^2);
        
        if abs(P1-P2)<1e-9
            D1=d/2;            
        else
            D1=[(-P1*d+sqrt(P1^3*P2+P2^3*P1-2*P1^2*P2^2+d^2*P1*P2))/(P2-P1),...
                (-P1*d-sqrt(P1^3*P2+P2^3*P1-2*P1^2*P2^2+d^2*P1*P2))/(P2-P1),...
                (+P1*d+sqrt(P1^3*P2+P2^3*P1-2*P1^2*P2^2+d^2*P1*P2))/(P2-P1),...
                (+P1*d-sqrt(P1^3*P2+P2^3*P1-2*P1^2*P2^2+d^2*P1*P2))/(P2-P1)];  
            
            D1=D1(D1>0);
            if isempty(D1)
                return;
            end
        end
        
        xi=[]; yi=[]; dd=[]; dw=[];
        for i=1:length(D1)
            
            if sign(sin(ksi1-AD1))==1
                XL1=xs1-D1(i)*cos(AD1); YL1=ys1-D1(i)*sin(AD1);
                XR1=xs1+D1(i)*cos(AD1); YR1=ys1+D1(i)*sin(AD1);                
            elseif sign(sin(ksi1-AD1))==-1
                XL1=xs1+D1(i)*cos(AD1); YL1=ys1+D1(i)*sin(AD1);
                XR1=xs1-D1(i)*cos(AD1); YR1=ys1-D1(i)*sin(AD1);   
            end
            
            M=[[-sin(ac), cos(ac)];...
               [-sin(AD1-pi/2), cos(AD1-pi/2)]];
            m=[ -sin(ac)*xq + cos(ac)*yq;...
                -sin(AD1-pi/2)*XL1 + cos(AD1-pi/2)*YL1];

            DM=det(M);
            if abs(DM)<1e-12
                return;
            end
            temp=M\m;
            xl1=temp(1);
            yl1=temp(2);
            
            M=[[-sin(ac), cos(ac)];...
               [-sin(AD1-pi/2), cos(AD1-pi/2)]];
            m=[ -sin(ac)*xq + cos(ac)*yq;...
                -sin(AD1-pi/2)*XR1 + cos(AD1-pi/2)*YR1];

            DM=det(M);
            if abs(DM)<1e-12
                return;
            end
            temp=M\m;
            xr1=temp(1);
            yr1=temp(2);
            
            DL1=sqrt((XL1-xl1)^2+(YL1-yl1)^2);
            DR1=sqrt((XR1-xr1)^2+(YR1-yr1)^2);            
            DFL=sqrt((XF1-xl1)^2+(YF1-yl1)^2);
            DFR=sqrt((XF1-xr1)^2+(YF1-yr1)^2);
            
            if abs(DFL-DL1)<1e-9 && abs(DFR-DR1)>1e-9
                xc=xl1; yc=yl1; dc1=-D1(i);
                stfd=1;
            elseif abs(DFL-DL1)>1e-9 && abs(DFR-DR1)<1e-9
                xc=xr1; yc=yr1; dc1=+D1(i);
                stfd=1;
            else
                stfd=0;
            end

            if stfd
                M=[[-sin(AD2-pi/2) , cos(AD2-pi/2)  ];...
                   [-sin(AD2) , cos(AD2) ]];
                m=[ -sin(AD2-pi/2)*XF2 + cos(AD2-pi/2)*YF2  ;...
                    -sin(AD2)*xc + cos(AD2)*yc ];

                DM=det(M);
                if abs(DM)<1e-12
                    return;
                end
                temp=M\m;
                xt=temp(1);
                yt=temp(2);

                sc2=sign((xc-xt)*sin(ksi2) - (yc-yt)*cos(ksi2));
                dc2=sc2*sqrt((xc-xt)^2+(yc-yt)^2);                 
                
                if cos(ksi1-a1)<0 && dp1>0
                    if dc1<dp1
                        status1=1;
                    else
                        status1=0;
                    end
                elseif cos(ksi1-a1)>0 && dp1>0
                    if dc1>dp1
                        status1=1;
                    else
                        status1=0;
                    end
                elseif cos(ksi1-a1)<0 && dp1<0    
                    if dc1>dp1
                        status1=1;
                    else
                        status1=0;
                    end
                elseif cos(ksi1-a1)>0 && dp1<0  
                    if dc1<dp1
                        status1=1;
                    else
                        status1=0;
                    end
                end
                
                if cos(ksi2-a2)<0 && dp2>0
                    if dc2<dp2
                        status2=1;
                    else
                        status2=0;
                    end
                elseif cos(ksi2-a2)>0 && dp2>0
                    if dc2>dp2
                        status2=1;
                    else
                        status2=0;
                    end
                elseif cos(ksi2-a2)<0 && dp2<0    
                    if dc2>dp2
                        status2=1;
                    else
                        status2=0;
                    end
                elseif cos(ksi2-a2)>0 && dp2<0  
                    if dc2<dp2
                        status2=1;
                    else
                        status2=0;
                    end
                end
                
                if status1==1 && status2==1                    
                    xi=[xi xc];
                    yi=[yi yc];
                    dd=[dd abs(dc1-dp1)];
                    dw=[dw abs(dc2-dp2)];            
                end
            end            
        end
        if isempty(dd) 
            return;
        else
            [dc,it]=min(dd);
            dm=dw(it);
            xc=xi(it);
            yc=yi(it);        
            
            am=atan2(0.5*sin(a1)+0.5*sin(a2),0.5*cos(a1)+0.5*cos(a2));
            if sign(cos(ac-am))==-1
                ac=ac-pi;
            end
            if ac<0
                ac=ac+2*pi;
            end

            tc=2;
            pc=[xq yq xq yq K1 K2];  
            status = 1;
        end
    end    
    
                
return;

