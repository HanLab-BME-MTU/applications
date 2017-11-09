function [dc,xc,yc,ac,tc,pc,status] = LineParabLine(xo1,yo1,ao1,to1,prm1,xo2,yo2,ao2,to2,prm2,XB,YB,AL,AR)

    pi=3.1415926535897932384626433832795;
    status=0; dc=[]; xc=[]; yc=[]; ac=[]; tc=[]; pc=[];

    if to1==0
        ksi=prm1(3); K=prm1(5); N=prm1(6);
        P=prm2(5);
        AN=ao2;
        xp=xo1; yp=yo1; %param=prm1;
        aa=ao1;        
    else 
        ksi=prm2(3); K=prm2(5); N=prm2(6);
        P=prm1(5);
        AN=ao1;
        xp=xo2; yp=yo2; %param=prm2;
        aa=ao2;
    end
            
    XF=XB(K); YF=YB(K);
    XD=XB(N); YD=YB(N); 
    XI=XB(P); YI=YB(P);
    
    if abs(cos(ksi-AL(N)))<1e-9 && abs(cos(ksi-AR(N)))>1e-9
        AD=AL(N);
    elseif abs(cos(ksi-AL(N)))>1e-9 && abs(cos(ksi-AR(N)))<1e-9
        AD=AR(N);
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else
        return;
    end
    
    if abs(cos(AN-AL(P)))<1e-9 && abs(cos(AN-AR(P)))>1e-9
        AI=AL(P);
    elseif abs(cos(AN-AL(P)))>1e-9 && abs(cos(AN-AR(P)))<1e-9
        AI=AR(P);
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else
        return;
    end
    
    %{
    M=[[-sin(AD), cos(AD)];...
       [-sin(AD-pi/2), cos(AD-pi/2)]];
    m=[ -sin(AD)*XD + cos(AD)*YD;...
        -sin(AD-pi/2)*XF + cos(AD-pi/2)*YF];

    DM=det(M);
    if abs(DM)<1e-12
        return;
    end
    temp=M\m;
    xs=temp(1);
    ys=temp(2);
    xo=(XF+xs)/2;
    yo=(YF+ys)/2;
    %}
  
    if K==P  
        xq=0; yq=0;
        tc=1;
        
        ac=atan2(sin(AD)+sin(AI),cos(AD)+cos(AI));
        if ac<0
           ac=ac+2*pi;
        end
                
        if sign(sin(ac-AI))==1 && sign(sin(ac-AD))==-1
            N3=N; M3=K;
        elseif sign(sin(ac-AI))==-1 && sign(sin(ac-AD))==1
            N3=K; M3=N;
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else
            return;       
        end        
    
        M=[[-sin(AI) , cos(AI)  ];...
           [-sin(AD) , cos(AD) ]];
        m=[ -sin(AI)*XF + cos(AI)*YF  ;...
            -sin(AD)*XD + cos(AD)*YD ];
        
        DM=det(M);
        if abs(DM)<1e-12
            return;
        end
        temp=M\m;
        xm=temp(1);
        ym=temp(2);
        
        M=[[-sin(AI-pi/2) , cos(AI-pi/2)  ];...
           [-sin(ac) , cos(ac) ]];
        m=[ -sin(AI-pi/2)*XF + cos(AI-pi/2)*YF  ;...
            -sin(ac)*xm + cos(ac)*ym ];
        
        DM=det(M);
        if abs(DM)<1e-12
            return;
        end
        temp=M\m;
        xc=temp(1);
        yc=temp(2);
                
    else
        tc=2;
        xq=(XI+XF)/2; yq=(YI+YF)/2;
        
        psi=atan2(YI-YF,XI-XF);
        if psi<0
           psi=psi+2*pi;
        end
        
        ac=psi-pi/2;
        if ac<0
            ac=ac+2*pi;              
        end    
        
        N3=K; M3=P;    
        
        M=[[-sin(psi-pi/2) , cos(psi-pi/2)  ];...
           [-sin(AD) , cos(AD) ]];
        m=[ -sin(psi-pi/2)*xq + cos(psi-pi/2)*yq  ;...
            -sin(AD)*XD + cos(AD)*YD ];
        
        DM=det(M);
        if abs(DM)<1e-12           
            return;
        end
        temp=M\m;
        xm=temp(1);
        ym=temp(2);
    
        M=[[-sin(AI-pi/2) , cos(AI-pi/2)  ];...
           [-sin(psi-pi/2) , cos(psi-pi/2) ]];
        m=[ -sin(AI-pi/2)*XI + cos(AI-pi/2)*YI  ;...
            -sin(psi-pi/2)*xq + cos(psi-pi/2)*yq ];
        
        DM=det(M);
        if abs(DM)<1e-12
            return;
        end
        temp=M\m;
        xc=temp(1);
        yc=temp(2);               
    end    
    
    M=[[-sin(AD-pi/2) , cos(AD-pi/2)  ];...
       [-sin(AD) , cos(AD) ]];
    m=[ -sin(AD-pi/2)*XF + cos(AD-pi/2)*YF  ;...
        -sin(AD)*xc + cos(AD)*yc ];

    DM=det(M);
    if abs(DM)<1e-12
        return;
    end
    temp=M\m;
    xt=temp(1);
    yt=temp(2);

    sc=sign((xc-xt)*sin(ksi) - (yc-yt)*cos(ksi));
    dc=sc*sqrt((xc-xt)^2+(yc-yt)^2);        

    M=[[-sin(AD-pi/2) , cos(AD-pi/2)  ];...
       [-sin(AD) , cos(AD) ]];
    m=[ -sin(AD-pi/2)*XF + cos(AD-pi/2)*YF  ;...
        -sin(AD)*xp + cos(AD)*yp ];

    DM=det(M);
    if abs(DM)<1e-12
        return;
    end
    temp=M\m;
    xt=temp(1);
    yt=temp(2);

    sp=sign((xp-xt)*sin(ksi) - (yp-yt)*cos(ksi));
    dp=sp*sqrt((xp-xt)^2+(yp-yt)^2); 

    if cos(ksi-aa)<0 && dp>0
        if dc<dp
            status=1;
        end
    elseif cos(ksi-aa)>0 && dp>0
        if dc>dp
            status=1;
        end
    elseif cos(ksi-aa)<0 && dp<0    
        if dc>dp
            status=1;
        end
    elseif cos(ksi-aa)>0 && dp<0  
        if dc<dp
            status=1;
        end
    end

    dc=abs(dc-dp); 
    pc=[xm ym xq yq N3 M3];      
    
    if cos(ksi-pi/2-ac)*cos(ksi-pi/2-aa)<0     
        ac=ac-pi;
    end
    if ac<0
        ac=ac+2*pi;
    end
    
return;  
