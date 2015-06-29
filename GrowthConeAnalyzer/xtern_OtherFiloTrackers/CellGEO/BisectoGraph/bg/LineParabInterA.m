function [X,Y,D,stfd]=LineParabInterA(XF1,YF1,PH1,XR1,YR1,PH2,XR2,YR2)

    CPH1=cos(PH1);
    SPH1=sin(PH1);
    
    XP1=XR1+(XF1-XR1)*CPH1^2+(YF1-YR1)*CPH1*SPH1;
    YP1=YR1+(XF1-XR1)*CPH1*SPH1+(YF1-YR1)*SPH1^2;

    XO1=(XF1+XP1)/2;
    YO1=(YF1+YP1)/2;
    
    xr=(XR2-XO1)*CPH1+(YR2-YO1)*SPH1;
    yr=(YR2-YO1)*CPH1-(XR2-XO1)*SPH1;
    yf1=(YF1-YO1)*CPH1-(XF1-XO1)*SPH1;
    phi=PH2-PH1;
    P1=2*yf1;
        
    cphi=cos(phi);
    sphi=sin(phi);
    
    A=cphi/(2*P1);
    B=-sphi;
    C=xr*sphi-yr*cphi;
    
    xs=roots([A,B,C]); 
    x=zeros(2,1); y=x; d=x;

    n=0;
    for i=1:2
        if imag(xs(i))==0
            n=n+1;
            x(n)=XO1+xs(i)*CPH1-(xs(i)^2/(2*P1))*SPH1;
            y(n)=YO1+(xs(i)^2/(2*P1))*CPH1+xs(i)*SPH1;
            if P1>=0
                d(n)=xs(i);
            else
                d(n)=-xs(i);
            end
            %f(n)=fparr(x(n),y(n),XF1,YF1,XR1,YR1,CPH1,SPH1);            
        end
    end
    X=x(1:n); Y=y(1:n); D=d(1:n); 
    stfd=n>0;
    
        
        
    
    

return;



    M=[[-sin(ADo), cos(ADo)];...
       [-sin(ADo-pi/2), cos(ADo-pi/2)]];
    m=[ -sin(ADo)*XDo + cos(ADo)*YDo;...
        -sin(ADo-pi/2)*XF + cos(ADo-pi/2)*YF];

    %DM=det(M);
    temp=M\m;
    xs=temp(1);
    ys=temp(2);
    xo=(XF+xs)/2;
    yo=(YF+ys)/2;
    P=sqrt((XF-xs)^2+(YF-ys)^2);
        
    M=[[-sin(aco), cos(aco)];...
       [-sin(ADo-pi/2), cos(ADo-pi/2)]];
    m=[ -sin(aco)*xco + cos(aco)*yco;...
        -sin(ADo-pi/2)*XF + cos(ADo-pi/2)*YF];

    %DM=det(M);
    temp=M\m;
    xu=temp(1);
    yu=temp(2);
    
    M=[[-sin(aco), cos(aco)];...
       [-sin(ADo), cos(ADo)]];
    m=[ -sin(aco)*xco + cos(aco)*yco;...
        -sin(ADo)*xo + cos(ADo)*yo];

    %DM=det(M);
    temp=M\m;
    xv=temp(1);
    yv=temp(2);
    
    xn=xo+(P/2)*cos(ADo);
    yn=yo+(P/2)*sin(ADo);    
    if double((xn-xo)*(YF-yo)-(yn-yo)*(XF-xo)) < 0
        ADo=pi+ADo;
        xn=xo+(P/2)*cos(ADo);
        yn=yo+(P/2)*sin(ADo);        
    end
    
    A=sqrt((xu-xo)^2+(yu-yo)^2);
    B=sqrt((xv-xo)^2+(yv-yo)^2);
    
    if double((xn-xo)*(yu-yo)-(yn-yo)*(xu-xo) ) > 0
        stfd=1;
        D1 = - P*A/B + sqrt (P*P*A*A/(B*B) + 2*P*A);
        D2 = - P*A/B - sqrt (P*P*A*A/(B*B) + 2*P*A);
        
        if double((XF-xo)*(yv-yo)-(YF-yo)*(xv-xo) ) > 0
            D1=-D1;
            D2=-D2;
        end
        
        X1=xo+D1*cos(ADo);
        Y1=yo+D1*sin(ADo);
        X2=xo+D2*cos(ADo);
        Y2=yo+D2*sin(ADo);
    else
        
        if double(P*P*A*A/(B*B)-2*P*A)<0
            X=[]; Y=[]; D=[]; stfd=0;            
        else
            stfd=1;
            D1 = P*A/B + sqrt (P*P*A*A/(B*B) - 2*P*A);
            D2 = P*A/B - sqrt (P*P*A*A/(B*B) - 2*P*A);

            if double((XF-xo)*(yv-yo)-(YF-yo)*(xv-xo) ) > 0
                D1=-D1;
                D2=-D2;
            end

            X1=xo+D1*cos(ADo);
            Y1=yo+D1*sin(ADo);
            X2=xo+D2*cos(ADo);
            Y2=yo+D2*sin(ADo);
        end
    end
    
    if stfd
    
        M=[[-sin(aco), cos(aco)];...
           [-sin(ADo-pi/2), cos(ADo-pi/2)]];
        m=[ -sin(aco)*xco + cos(aco)*yco;...
            -sin(ADo-pi/2)*X1 + cos(ADo-pi/2)*Y1];

        %DM=det(M);
        temp=M\m;
        xc1=temp(1);
        yc1=temp(2);

        M=[[-sin(aco), cos(aco)];...
           [-sin(ADo-pi/2), cos(ADo-pi/2)]];
        m=[ -sin(aco)*xco + cos(aco)*yco;...
            -sin(ADo-pi/2)*X2 + cos(ADo-pi/2)*Y2];

        %DM=det(M);
        temp=M\m;
        xc2=temp(1);
        yc2=temp(2);

        X=[double(xc1) double(xc2)];
        Y=[double(yc1) double(yc2)];
        D=[double(D1)  double(D2)];
             
    end
    
return;
    