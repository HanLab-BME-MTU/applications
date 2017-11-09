function param = LineLineParab(x1,y1,a1,t1,prm1,x2,y2,a2,t2,prm2,XB,YB,AL,AR)

    pi=3.1415926535897932384626433832795;
    
    param = [];
         
    if t1==1 
        AQ=a1; AP=a2;
        N1=prm1(5); N2=prm1(6);
        K=prm2(5);
    else 
        AQ=a2; AP=a1;
        N1=prm2(5); N2=prm2(6);
        K=prm1(5);
    end
    
    if N1==N2        
        AAL=AL(N1);
        AAR=AR(N1);
    else
        AA1=AL(N1);
        AA2=AR(N1);
        AA3=AL(N2);
        AA4=AR(N2);

        AM13=atan2(sin(AA1)+sin(AA3),cos(AA1)+cos(AA3));
        if AM13<0
            AM13=AM13+2*pi;
        end
        AM14=atan2(sin(AA1)+sin(AA4),cos(AA1)+cos(AA4));
        if AM14<0
            AM14=AM14+2*pi;
        end
        AM23=atan2(sin(AA2)+sin(AA3),cos(AA2)+cos(AA3));
        if AM23<0
            AM23=AM23+2*pi;
        end
        AM24=atan2(sin(AA2)+sin(AA4),cos(AA2)+cos(AA4));
        if AM24<0
            AM24=AM24+2*pi;
        end

        if abs(sin(AQ-AM13))<1e-9
            AAL=AA1; AAR=AA3;
        elseif abs(sin(AQ-AM14))<1e-9
            AAL=AA1; AAR=AA4;
        elseif abs(sin(AQ-AM23))<1e-9
            AAL=AA2; AAR=AA3;
        elseif abs(sin(AQ-AM24))<1e-9
            AAL=AA2; AAR=AA4;
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else
            return;
        end    
    end
            
    if abs(cos(AAL-AP))<1e-9
        N=N2;
        AD=AAR;
    elseif abs(cos(AAR-AP))<1e-9
        N=N1;
        AD=AAL;
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else
        return;    
    end        
      
    XD=XB(N); YD=YB(N);
    XF=XB(K); YF=YB(K);
    
    M=[[-sin(AD)     , cos(AD)      ];...
       [-sin(AD-pi/2), cos(AD-pi/2) ]];
    m=[ -sin(AD)*XD + cos(AD)*YD  ;...
        -sin(AD-pi/2)*XF + cos(AD-pi/2)*YF];

    DM=det(M);
    if abs(DM)<1e-12
        return;
    end
    temp=M\m;
    XS=temp(1);
    YS=temp(2);

    XO=(XF+XS)/2;
    YO=(YF+YS)/2;
    PO=sqrt((XF-XS)^2+(YF-YS)^2);
    AO=atan2(YF-YS,XF-XS);
    if AO<0
        AO=AO+2*pi;
    end
    
    ksi=AD-pi/2;
    if ksi<0
        ksi=ksi+2*pi;
    end
    
    if cos(ksi-AO)<0
        ksi=AD+pi/2;
        if ksi>(2*pi);
            ksi=ksi-2*pi;
        end
    end
        
    param=[XO YO ksi PO K N];
    
return;