function status = ParabPlotX(x1,y1,prm1,x2,y2,XB,YB,AL,AR,col)

    status = 0;
    ksi=prm1(3);
    K=prm1(5); N=prm1(6);
    
    XF=XB(K); YF=YB(K);
    XD=XB(N); YD=YB(N);
    
    if abs(cos(ksi-AL(N)))<1e-9 && abs(cos(ksi-AR(N)))>1e-9
        AD=AL(N);
    elseif abs(cos(ksi-AL(N)))>1e-9 && abs(cos(ksi-AR(N)))<1e-9
        AD=AR(N);    
    end
    
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
    p=sqrt((XF-xs)^2+(YF-ys)^2);
    
    M=[[-sin(AD-pi/2), cos(AD-pi/2)      ];...
       [-sin(AD) ,     cos(AD) ]];
    m=[ -sin(AD-pi/2)*x1 + cos(AD-pi/2)*y1 ;...
        -sin(AD)*xo + cos(AD)*yo];

    DM=det(M);
    if abs(DM)<1e-12
        return;
    end
    temp=M\m;
    xp1=temp(1);
    yp1=temp(2);
    
    M=[[-sin(AD-pi/2), cos(AD-pi/2)      ];...
       [-sin(AD),      cos(AD) ]];
    m=[ -sin(AD-pi/2)*x2 + cos(AD-pi/2)*y2 ;...
        -sin(AD)*xo + cos(AD)*yo];

    DM=det(M);
    if abs(DM)<1e-12
        return;
    end
    temp=M\m;
    xp2=temp(1);
    yp2=temp(2);
         
    N=20;
    xr=linspace(xp1,xp2,N);
    yr=linspace(yp1,yp2,N);
    dr=((xr-xo).^2+(yr-yo).^2)/(2*p);
    X=xr+dr*cos(ksi);
    Y=yr+dr*sin(ksi);
        
    %plot(X,Y,'g','LineWidth',2);
    plot(X,Y,'Color',col,'LineWidth',1);
    
      