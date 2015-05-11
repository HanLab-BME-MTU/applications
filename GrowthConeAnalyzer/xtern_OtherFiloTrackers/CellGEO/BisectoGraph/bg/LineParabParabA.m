function [dc,xc,yc,ac,tc,pc,status] = LineParabParabA(X,Y,phi,prm,xa,ya,aa,param,XB,YB,AL,AR)

    dc=[]; xc=[]; yc=[]; ac=[]; tc=[]; pc=[]; status=0;
    
    ksi=param(3); 

    K=param(5); %N=param(6);
    N1=prm(5); M1=prm(6);
    
    al=AL(N1);
    ar=AR(M1);
    
    if abs(cos(ksi-al))<1e-9 && abs(cos(ksi-ar))>1e-9
        ADo=al; ADn=ar;
        No=N1; Nn=M1;        
    elseif abs(cos(ksi-al))>1e-9 && abs(cos(ksi-ar))<1e-9
        ADo=ar; ADn=al;
        No=M1; Nn=N1;           
    else
        al=AR(N1);
        ar=AL(M1);
    
        if abs(cos(ksi-al))<1e-9 && abs(cos(ksi-ar))>1e-9
            ADo=al; ADn=ar;
            No=N1; Nn=M1;        
        elseif abs(cos(ksi-al))>1e-9 && abs(cos(ksi-ar))<1e-9
            ADo=ar; ADn=al;
            No=M1; Nn=N1;
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else
            return;        
        end
    end
    
    XDo=XB(No); YDo=YB(No);
    XDn=XB(Nn); YDn=YB(Nn);   
    XF=XB(K); YF=YB(K);
     
    
    
    
    
    M=[[-sin(ADo), cos(ADo)];...
       [-sin(ADo-pi/2), cos(ADo-pi/2)]];
    m=[ -sin(ADo)*XDo + cos(ADo)*YDo;...
        -sin(ADo-pi/2)*XF + cos(ADo-pi/2)*YF];

    DM=det(M);
    if abs(DM)<1e-12
        return;
    end
    temp=M\m;
    xs=temp(1);
    ys=temp(2);
    xo=(XF+xs)/2;
    yo=(YF+ys)/2;
    %P=sqrt((XF-xs)^2+(YF-ys)^2);
    
    
    %---------------
    
    aco=atan2(sin(ADo)+sin(ADn),cos(ADo)+cos(ADn));
    if aco<0
        aco=aco+2*pi;
    end
    
    M=[[-sin(aco-pi/2) , cos(aco-pi/2)  ];...
       [-sin(ADn) , cos(ADn) ]];
    m=[ -sin(aco-pi/2)*XDo + cos(aco-pi/2)*YDo  ;...
        -sin(ADn)*XDn + cos(ADn)*YDn ];

    DM=det(M);
    if abs(DM)<1e-12
        return;
    end
    temp=M\m;
    xco=(XDo+temp(1))/2;
    yco=(YDo+temp(2))/2;
    
    
    M=[[-sin(ADo-pi/2) , cos(ADo-pi/2)  ];...
       [-sin(ADo) , cos(ADo) ]];
    m=[ -sin(ADo-pi/2)*XF + cos(ADo-pi/2)*YF  ;...
        -sin(ADo)*xa + cos(ADo)*ya ];

    DM=det(M);
    if abs(DM)<1e-12
        return;
    end
    temp=M\m;
    xt=temp(1);
    yt=temp(2);

    sp=sign((xa-xt)*sin(ksi) - (ya-yt)*cos(ksi));
    dp=sp*sqrt((xa-xt)^2+(ya-yt)^2); 

    [XC,YC,DC,stfd]=LineParabInterA(XF,YF,ADo,XDo,YDo,aco,xco,yco);
    
    xi=[]; yi=[]; dd=[];
    for i=1:2
        if stfd
            xc=XC(i); yc=YC(i); dc=DC(i);
            if cos(ksi-aa)<0 && dp>0
                if dc<dp
                    xi=[xi xc];
                    yi=[yi yc];
                    dd=[dd abs(dc-dp)];
                end
            elseif cos(ksi-aa)>0 && dp>0
                if dc>dp
                    xi=[xi xc];
                    yi=[yi yc];
                    dd=[dd abs(dc-dp)];
                end
            elseif cos(ksi-aa)<0 && dp<0    
                if dc>dp
                    xi=[xi xc];
                    yi=[yi yc];
                    dd=[dd abs(dc-dp)];
                end
            elseif cos(ksi-aa)>0 && dp<0  
                if dc<dp
                    xi=[xi xc];
                    yi=[yi yc];
                    dd=[dd abs(dc-dp)];
                end
            end
        end
    end

    if isempty(dd) 
        return;
    else
        [dc,it]=min(dd);
        xc=xi(it);
        yc=yi(it);    
        tc=0;
        
        M=[[-sin(ADn), cos(ADn)];...
           [-sin(ADn-pi/2), cos(ADn-pi/2)]];
        m=[ -sin(ADn)*XDn + cos(ADn)*YDn;...
            -sin(ADn-pi/2)*XF + cos(ADn-pi/2)*YF];

        DM=det(M);
        if abs(DM)<1e-12
            return;
        end
        temp=M\m;
        nxs=temp(1);
        nys=temp(2);        
        nxo=(XF+nxs)/2;
        nyo=(YF+nys)/2;
        np=sqrt((XF-nxs)^2+(YF-nys)^2);
        nksi=atan2(YF-nys,XF-nxs);
        if nksi<0
            nksi=nksi+2*pi;
        end
        
        %----------------------------------------------------------------------
      
        M=[[-sin(ADo) ,      cos(ADo)  ];...
           [-sin(ADo-pi/2) , cos(ADo-pi/2) ]];
        m=[ -sin(ADo)*xo + cos(ADo)*yo  ;...
            -sin(ADo-pi/2)*xc + cos(ADo-pi/2)*yc ];

        DM=det(M);
        if abs(DM)<1e-12
            return;
        end
        temp=M\m;
        xp=temp(1);
        yp=temp(2);
        mxm=(xp+xo)/2;
        mym=(yp+yo)/2;
        
        aan=atan2(yc-mym,xc-mxm);
        if aan<0
            aan=aan+2*pi;
        end
        
        if sign(cos(aan-aa))==-1
            aan=aan-pi;
            if aan<0
                aan=aan+2*pi;
            end
        end
      
        %----------------------------------------------------------------------
  
        M=[[-sin(ADn) ,      cos(ADn)  ];...
           [-sin(ADn-pi/2) , cos(ADn-pi/2) ]];
        m=[ -sin(ADn)*nxo + cos(ADn)*nyo  ;...
            -sin(ADn-pi/2)*xc + cos(ADn-pi/2)*yc ];

        DM=det(M);
        if abs(DM)<1e-12
            return;
        end
        temp=M\m;
        nxp=temp(1);
        nyp=temp(2);
        nxm=(nxp+nxo)/2;
        nym=(nyp+nyo)/2;

        ac=atan2(yc-nym,xc-nxm);
        if ac<0
            ac=ac+2*pi;
        end

        %----------------------------------------------------------------------

        if  cos(aan-ac)<0 || cos(phi-ac)<0    
            ac=ac-pi;    
        end
        if ac<0
            ac=ac+2*pi;
        end

        tc=0;
        pc=[nxo nyo nksi np K Nn];
        status=1;
    end
    
return;


