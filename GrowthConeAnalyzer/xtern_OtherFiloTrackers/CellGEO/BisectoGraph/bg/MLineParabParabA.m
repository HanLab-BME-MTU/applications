function [dc,xc,yc,ac,tc,pc,status] = MLineParabParabA(X,Y,phi,prm,xa,ya,aa,param,XB,YB,AL,AR)

    dc=[]; xc=[]; yc=[]; ac=[]; tc=[]; pc=[]; status=0;
    %xo=param(1); yo=param(2); p=param(4);
    ksi=param(3); 

    K=param(5); N=param(6);
    K1=prm(5);  K2=prm(6);
    
    al=AL(N);
    ar=AR(N);
    
    if abs(cos(ksi-al))<1e-9 && abs(cos(ksi-ar))>1e-9
        ADo=al;      
    elseif abs(cos(ksi-al))>1e-9 && abs(cos(ksi-ar))<1e-9
        ADo=ar;            
    else
        al=AR(N);
        ar=AL(N);
        if abs(cos(ksi-al))<1e-9 && abs(cos(ksi-ar))>1e-9
            ADo=al;      
        elseif abs(cos(ksi-al))>1e-9 && abs(cos(ksi-ar))<1e-9
            ADo=ar;
        else
            disp('MLineParabParab - problem: AD is not consistent with ksi (line 21)');
            %abrakadabra
            return;
        end
    end
    
    if K==K1
        Ko=K1; Kn=K2;
    elseif K==K2
        Ko=K2; Kn=K1;
    else
        disp('MLineParabParab - problem: K~=K1 and K~=K2 (line 30)');
        %abrakadabra
        return;
    end
        
    XDo=XB(N);  YDo=YB(N);
    XFo=XB(Ko); YFo=YB(Ko);   
    XFn=XB(Kn); YFn=YB(Kn);
    xq=(XFn+XFo)/2; yq=(YFn+YFo)/2; 
     
    M=[[-sin(ADo), cos(ADo)];...
       [-sin(ADo-pi/2), cos(ADo-pi/2)]];
    m=[ -sin(ADo)*XDo + cos(ADo)*YDo;...
        -sin(ADo-pi/2)*XFo + cos(ADo-pi/2)*YFo];

    DM=det(M);
    if abs(DM)<1e-12
        disp(' ');
        disp('det(M) is too small (ParabParabInter, line 41)');
        disp(' ');
        return;
    end
    temp=M\m;
    xs=temp(1);
    ys=temp(2);
    xo=(XFo+xs)/2;
    yo=(YFo+ys)/2;
        
    M=[[-sin(ADo), cos(ADo)];...
       [-sin(ADo-pi/2), cos(ADo-pi/2)]];
    m=[ -sin(ADo)*XDo + cos(ADo)*YDo;...
        -sin(ADo-pi/2)*XFn + cos(ADo-pi/2)*YFn];

    DM=det(M);
    if abs(DM)<1e-12
        disp(' ');
        disp('det(M) is too small (ParabParabInter, line 41)');
        disp(' ');
        return;
    end
    temp=M\m;
    nxs=temp(1);
    nys=temp(2);
    nxo=(XFn+nxs)/2;
    nyo=(YFn+nys)/2;
    
    %xqs=(xs+nxs)/2;
    %yqs=(ys+nys)/2;
    %P=sqrt((XFo-xs)^2+(YFo-ys)^2);
        
    aco=-pi/2+atan2(YFn-YFo,XFn-XFo);
    if aco<0
        aco=aco+2*pi;
    end
    
    
    %{    
    M=[[-sin(ADo), cos(ADo)];...
       [-sin(aco), cos(aco)]];
    m=[ -sin(ADo)*XDo + cos(ADo)*YDo;...
        -sin(aco)*xq + cos(aco)*yq];

    DM=det(M);
    if abs(DM)<1e-12
        disp(' ');
        disp('det(M) is too small (ParabParabInter, line 41)');
        disp(' ');
        return;
    end
    temp=inv(M)*m;
    xm=temp(1);
    ym=temp(2);
    
    A=sqrt((xm-xs)^2+(ym-ys)^2);    
     T=sqrt((xq-xqs)^2+(yq-yqs)^2)/sqrt((xm-xqs)^2+(ym-yqs)^2);
    if (P^2*(T^2-1)-2*A*P*T)>0
        D=[-T*P+sqrt(P^2*(T^2-1)+2*A*P*T),...
           -T*P-sqrt(P^2*(T^2-1)+2*A*P*T),...
           +T*P+sqrt(P^2*(T^2-1)+2*A*P*T),...
           +T*P-sqrt(P^2*(T^2-1)+2*A*P*T),...
           -T*P+sqrt(P^2*(T^2-1)-2*A*P*T),...
           -T*P-sqrt(P^2*(T^2-1)-2*A*P*T),...
           +T*P+sqrt(P^2*(T^2-1)-2*A*P*T),...
           +T*P-sqrt(P^2*(T^2-1)-2*A*P*T)];
    else
        D=[-T*P+sqrt(P^2*(T^2-1)+2*A*P*T),...
           -T*P-sqrt(P^2*(T^2-1)+2*A*P*T),...
           +T*P+sqrt(P^2*(T^2-1)+2*A*P*T),...
           +T*P-sqrt(P^2*(T^2-1)+2*A*P*T)];
    end
    
    
    D=D(D>0);
    if isempty(D)
        disp('MLineParabParba - problem: no intersections found (line 127)');
        abra
    end   
    %----
    %}
    
    M=[[-sin(ADo-pi/2) , cos(ADo-pi/2)  ];...
       [-sin(ADo) , cos(ADo) ]];
    m=[ -sin(ADo-pi/2)*XFo + cos(ADo-pi/2)*YFo  ;...
        -sin(ADo)*xa + cos(ADo)*ya ];

    DM=det(M);
    if abs(DM)<1e-12
        disp(' ');
        disp('det(M) is too small (LineParabInter, line 61)');        
        disp(' ');
        return;
    end
    temp=M\m;
    xt=temp(1);
    yt=temp(2);

    sp=sign((xa-xt)*sin(ksi) - (ya-yt)*cos(ksi));
    dp=sp*sqrt((xa-xt)^2+(ya-yt)^2); 
    
    [XC,YC,DC,stfd]=LineParabInterA(XFo,YFo,ADo,XDo,YDo,aco,xq,yq);

    xi=[]; yi=[]; dd=[];
    for i=1:2%length(D)
        %{
        if sign(sin(ksi-ADo))==1
            XL=xs-D(i)*cos(ADo); YL=ys-D(i)*sin(ADo);
            XR=xs+D(i)*cos(ADo); YR=ys+D(i)*sin(ADo);                
        elseif sign(sin(ksi-ADo))==-1
            XL=xs+D(i)*cos(ADo); YL=ys+D(i)*sin(ADo);
            XR=xs-D(i)*cos(ADo); YR=ys-D(i)*sin(ADo);   
        end

        M=[[-sin(aco), cos(aco)];...
           [-sin(ADo-pi/2), cos(ADo-pi/2)]];
        m=[ -sin(aco)*xm + cos(aco)*ym;...
            -sin(ADo-pi/2)*XL + cos(ADo-pi/2)*YL];

        DM=det(M);
        if abs(DM)<1e-12
            disp(' ');
            disp('det(M) is too small (MParabParabLine, line 156)');
            disp(' ');
            return;
        end
        temp=inv(M)*m;
        xl=temp(1);
        yl=temp(2);

        M=[[-sin(aco), cos(aco)];...
           [-sin(ADo-pi/2), cos(ADo-pi/2)]];
        m=[ -sin(aco)*xm + cos(aco)*ym;...
            -sin(ADo-pi/2)*XR + cos(ADo-pi/2)*YR];

        DM=det(M);
        if abs(DM)<1e-12
            disp(' ');
            disp('det(M) is too small (MParabParabLine, line 172)');
            disp(' ');
            return;
        end
        temp=inv(M)*m;
        xr=temp(1);
        yr=temp(2);

        DL=sqrt((XL-xl)^2+(YL-yl)^2);
        DR=sqrt((XR-xr)^2+(YR-yr)^2);            
        DFL=sqrt((XFo-xl)^2+(YFo-yl)^2);
        DFR=sqrt((XFo-xr)^2+(YFo-yr)^2);

        if     abs(DFL-DL)<1e-9 && abs(DFR-DR)>1e-9
            xc=xl; yc=yl; dc=-D(i);   
            stfd=1;
        elseif abs(DFL-DL)>1e-9 && abs(DFR-DR)<1e-9
            xc=xr; yc=yr; dc=+D(i);               
            stfd=1;
        else
            stfd=0;
            %disp('ParabParabLine - problem: no intersections found (line 192)');  
            %return;
        end
        %}

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
        %disp('no the one');
        %[K1 N1 K2 N2]
        return;
    else
        %disp('it is the one');
        %[K1 N1 K2 N2]
        [dc,it]=min(dd);
        xc=xi(it);
        yc=yi(it);    
        
        %plot(xc,yc,'b*');
        %plot(xc,yc,'bo');
        
        tc=0;
        
        np=sqrt((XFn-nxs)^2+(YFn-nys)^2);
        nksi=atan2(YFn-nys,XFn-nxs);
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
            disp(' ');
            disp('det(M) is too small (LineParabInter, line 153)');        
            disp(' ');
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
      
        %!!!!!!!!!
        if sign(cos(aan-aa))==-1
            aan=aan-pi;
            if aan<0
                aan=aan+2*pi;
            end
        end
        %!!!!!!!!
        
        %----------------------------------------------------------------------
  
        M=[[-sin(ADo) ,      cos(ADo)  ];...
           [-sin(ADo-pi/2) , cos(ADo-pi/2) ]];
        m=[ -sin(ADo)*nxo + cos(ADo)*nyo  ;...
            -sin(ADo-pi/2)*xc + cos(ADo-pi/2)*yc ];

        DM=det(M);
        if abs(DM)<1e-12
            disp(' ');
            disp('det(M) is too small (LineParabInter, line 153)');        
            disp(' ');
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

        am=atan2(sin(phi)+sin(aan),cos(phi)+cos(aan));
        if am<0
           am=am+2*pi;
        end
        
        %if  cos(aan-ac)<0 || cos(phi-ac)<0    
        if sign(cos(am-ac))==-1
            ac=ac-pi;    
        end
        if ac<0
            ac=ac+2*pi;
        end
        
        %!!!!!!!!!!!!!!!
        %{
        if K==58 && N==279 && K1==56 && K2==58
            plot(xc,yc,'kp');
            line([xo xp],[yo yp],'Color','k');
            plot(mxm,mym,'ks');
            line([mxm xc],[mym yc],'Color','m');
            
            180*[phi aan am ac]/pi
        end
        %}
        %!!!!!!!!!!!!!!!!
        
        tc=0;
        pc=[nxo nyo nksi np Kn N];
        status=1;
    end
    
return;









