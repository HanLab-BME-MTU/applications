function [ac,pc] = MMLineLineLine(a1,prm1,a2,prm2,XB,YB,AL,AR)

    pi=3.1415926535897932384626433832795;
    N1=prm1(5); M1=prm1(6);
    N2=prm2(5); M2=prm2(6);
    
    XFL1=XB(N1); YFL1=YB(N1);    
    XFR1=XB(M1); YFR1=YB(M1); 
    XFL2=XB(N2); YFL2=YB(N2);
    XFR2=XB(M2); YFR2=YB(M2);
    
    if abs(XFL1-XFL2)+abs(YFL1-YFL2)<1e-9
        XFL=XFR1; YFL=YFR1; N3=M1;
        XFR=XFR2; YFR=YFR2; M3=M2;
    elseif abs(XFL1-XFR2)+abs(YFL1-YFR2)<1e-9
        XFL=XFR1; YFL=YFR1; N3=M1; 
        XFR=XFL2; YFR=YFL2; M3=N2;
    elseif abs(XFR1-XFL2)+abs(YFR1-YFL2)<1e-9
        XFL=XFL1; YFL=YFL1; N3=N1; 
        XFR=XFR2; YFR=YFR2; M3=M2; 
    elseif abs(XFR1-XFR2)+abs(YFR1-YFR2)<1e-9
        XFL=XFL1; YFL=YFL1; N3=N1; 
        XFR=XFL2; YFR=YFL2; M3=N2; 
    %------------------------ need to fix -----------------------------     
    else
        ac=[]; pc=[];
        disp('MMLineLineLine');
        return;
    %------------------------ need to fix -----------------------------     
    end
    
    xq3=(XFL+XFR)/2; yq3=(YFL+YFR)/2;   
    
    ac=-pi/2+atan2(YFL-YFR,XFL-XFR);
    if ac<0
        ac=ac+2*pi;
    end

    am=atan2(sin(a1)+sin(a2),cos(a1)+cos(a2));
    if am<0
       am=am+2*pi;
    end
    
    if cos(ac-am)<0
        ac=ac-pi;
    end    
    if ac<0
        ac=ac+2*pi;            
    end
    
    pc=[xq3 yq3 xq3 yq3 N3 M3];
    
return;