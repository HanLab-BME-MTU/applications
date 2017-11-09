function status = LinearityCheck(x1,y1,t1,prm1,x2,y2,XB,YB,MB)
    
   pi=3.1415926535897932384626433832795; 
   status=0;
   xq = prm1(3); yq = prm1(4);
   
   if abs(xq)+abs(yq)>1e-9
       if t1==0           
           K=prm1(5);
           xf=XB(K);
           yf=YB(K);
           st1=inpolygon(x1,y1,[XB'; XB(1)],[YB'; YB(1)]);
           zt1=DistToBoundX(x1,y1,XB,YB,MB);
       elseif t1==2
           ksi=atan2(y2-y1,x2-x1);
           
           K1=prm1(5); K2=prm1(6);
           XF1=XB(K1); YF1=YB(K1);
           XF2=XB(K2); YF2=YB(K2);
           
           xq=(XF1+XF2)/2;
           yq=(YF1+YF2)/2;
           
           st1=inpolygon(x1,y1,[XB'; XB(1)],[YB'; YB(1)]);
           zt1=DistToBoundX(x1,y1,XB,YB,MB);
           
           z2=sqrt((xq-x1).^2+(yq-y1).^2);           
           zq=sqrt(zt1^2-z2^2);
           
           xf=xq+zq*cos(ksi-pi/2);
           yf=yq+zq*sin(ksi-pi/2);
       end
       st2=inpolygon(x2,y2,[XB'; XB(1)],[YB'; YB(1)]);
       zt2=DistToBoundX(x2,y2,XB,YB,MB);
       
       if st1==0 || st2==0
           return;
       end
       
       zf1=sqrt((xf-x1).^2+(yf-y1).^2);
       zf2=sqrt((xf-x2).^2+(yf-y2).^2);
       if abs(zt1-zf1)+abs(zt2-zf2)<1e-6    %!!!!!!!!!!!!!!
          status=1;
       end
     
   elseif abs(xq)+abs(yq)<1e-9 && t1~=0   
       xm=(x1+x2)/2; ym=(y1+y2)/2;
       
       st1=inpolygon(x1,y1,[XB'; XB(1)],[YB'; YB(1)]);
       zt1=DistToBoundX(x1,y1,XB,YB,MB);
       stm=inpolygon(xm,ym,[XB'; XB(1)],[YB'; YB(1)]);
       ztm=DistToBoundX(xm,ym,XB,YB,MB);
       st2=inpolygon(x2,y2,[XB'; XB(1)],[YB'; YB(1)]);
       zt2=DistToBoundX(x2,y2,XB,YB,MB);

       if st1==0 || stm==0 || st2==0
           return;
       end

       ztc=(zt1+zt2)/2;
       if abs(ztm-ztc)<1e-9   %!!!!!!!!!!!!!!!!
           status=1;
       end  
   end

return;