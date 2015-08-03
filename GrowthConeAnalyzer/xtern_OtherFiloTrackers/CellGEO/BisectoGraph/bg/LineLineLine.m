function [ac,pc] = LineLineLine(a1,prm1,a2,prm2,XB,YB,AL,AR)
 
    pi=3.1415926535897932384626433832795;
    N1=prm1(5); M1=prm1(6);
    N2=prm2(5); M2=prm2(6);
    
    ac=[]; pc=[];
    
    %----------------------------------------------------------
    aNL=AL(N1);
    aNR=AR(N1);
    aML=AL(M1);
    aMR=AR(M1);
    
    aLL=atan2(sin(aNL)+sin(aML),cos(aNL)+cos(aML));
    if aLL<0
       aLL=aLL+2*pi;
    end
    aLR=atan2(sin(aNL)+sin(aMR),cos(aNL)+cos(aMR));
    if aLR<0
       aLR=aLR+2*pi;
    end
    aRL=atan2(sin(aNR)+sin(aML),cos(aNR)+cos(aML));
    if aRL<0
       aRL=aRL+2*pi;
    end
    aRR=atan2(sin(aNR)+sin(aMR),cos(aNR)+cos(aMR));
    if aRR<0
       aRR=aRR+2*pi;
    end
    
    [tv ti]=min(abs([a1 a1 a1 a1]-[aLL aLR aRL aRR]));
    
    if ti==1
        bL1=aNL; bR1=aML;
    elseif ti==2
        bL1=aNL; bR1=aMR;
    elseif ti==3
        bL1=aNR; bR1=aML;
    elseif ti==4
        bL1=aNR; bR1=aMR;
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else
        return;    
    end
    
    %----------------------------------------------------------
    aNL=AL(N2);
    aNR=AR(N2);
    aML=AL(M2);
    aMR=AR(M2);
    
    aLL=atan2(sin(aNL)+sin(aML),cos(aNL)+cos(aML));
    if aLL<0
       aLL=aLL+2*pi;
    end
    aLR=atan2(sin(aNL)+sin(aMR),cos(aNL)+cos(aMR));
    if aLR<0
       aLR=aLR+2*pi;
    end
    aRL=atan2(sin(aNR)+sin(aML),cos(aNR)+cos(aML));
    if aRL<0
       aRL=aRL+2*pi;
    end
    aRR=atan2(sin(aNR)+sin(aMR),cos(aNR)+cos(aMR));
    if aRR<0
       aRR=aRR+2*pi;
    end
    
    [tv ti]=min(abs([a2 a2 a2 a2]-[aLL aLR aRL aRR]));
    
    if ti==1
        bL2=aNL; bR2=aML;
    elseif ti==2
        bL2=aNL; bR2=aMR;
    elseif ti==3
        bL2=aNR; bR2=aML;
    elseif ti==4
        bL2=aNR; bR2=aMR;
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else
        return;    
    end
    
    %----------------------------------------------------------
    bLL=abs(sin(bL1-bL2));
    bLR=abs(sin(bL1-bR2));
    bRL=abs(sin(bR1-bL2));
    bRR=abs(sin(bR1-bR2));
    
    [tv ti]=min([bLL bLR bRL bRR]);
    
    if ti==1
        AL3=bR1; AR3=bR2;
        KL3=M1;  KR3=M2;
    elseif ti==2
        AL3=bR1; AR3=bL2;
        KL3=M1;  KR3=N2;
    elseif ti==3
        AL3=bL1; AR3=bR2;
        KL3=N1;  KR3=M2;
    elseif ti==4
        AL3=bL1; AR3=bL2;
        KL3=N1;  KR3=N2;
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else
        return;    
    end
    
    XL3=XB(KL3); YL3=YB(KL3);
    XR3=XB(KR3); YR3=YB(KR3);
    
    
%--------------------------------------------------------------------------
    
    M=[[-sin(AL3) , cos(AL3)      ];...
       [-sin(AR3) , cos(AR3) ]];
    m=[ -sin(AL3)*XL3 + cos(AL3)*YL3  ;...
        -sin(AR3)*XR3 + cos(AR3)*YR3];

    DM=det(M);
    if abs(DM)<1e-12
        return;
    end
    temp=M\m;
    xt=temp(1);
    yt=temp(2);
    
    
    ac=atan2(sin(AL3)+sin(AR3),cos(AL3)+cos(AR3));
    if ac<0
       ac=ac+2*pi;
    end
    
    pc=[xt yt 0 0 KL3 KR3];
