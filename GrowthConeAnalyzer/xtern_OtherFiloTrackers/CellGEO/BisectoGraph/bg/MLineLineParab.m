function param = MLineLineParab(x1,y1,a1,t1,prm1,x2,y2,a2,t2,prm2,XB,YB,AL,AR)

    pi=3.1415926535897932384626433832795;
    if t1==2 
        N=prm1(5); K=prm1(6);
        P=prm2(5);
        AN=a2;
    else 
        N=prm2(5); K=prm2(6);
        P=prm1(5);
        AN=a1;        
    end
    
    if K==P
        XF=XB(N); YF=YB(N);
        NF=N;
    elseif N==P
        XF=XB(K); YF=YB(K);
        NF=K;
    %------------------------ need to fix -----------------------------    
    else
        param = [];
        disp('MLineLineParab');
        return;        
    %------------------------ need to fix -----------------------------    
    end 
    
    XI=XB(P); YI=YB(P);
    
    if abs(cos(AN-AL(P)))<1e-9 && abs(cos(AN-AR(P)))>1e-9
        AI=AL(P);
    elseif abs(cos(AN-AL(P)))>1e-9 && abs(cos(AN-AR(P)))<1e-9
        AI=AR(P);
    end    
    
    M=[[-sin(AI)      , cos(AI)      ];...
       [-sin(AI-pi/2) , cos(AI-pi/2) ]];
    m=[ -sin(AI)*XI + cos(AI)*YI  ;...
        -sin(AI-pi/2)*XF + cos(AI-pi/2)*YF];

    DM=det(M);
    if abs(DM)<1e-12        
        return;
    end
    %temp=inv(M)*m;
    temp=M\m;
    XS=temp(1);
    YS=temp(2);

    XO=(XF+XS)/2;
    YO=(YF+YS)/2;

    p=sqrt((XF-XS)^2+(YF-YS)^2);
    param=[XO YO AN p NF P];

return;