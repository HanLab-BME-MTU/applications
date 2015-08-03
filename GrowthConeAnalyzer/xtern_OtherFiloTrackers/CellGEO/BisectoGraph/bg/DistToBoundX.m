function db=DistToBoundX(xo,yo,XB,YB,MB)

    if xo>MB(1) || xo<MB(2) || yo>MB(3) || yo<MB(4)
        db=0;
        disp('wtf');
        return;
    end
    
    LB=length(XB);
    e=zeros(1,LB);
    for k=1:(LB-1)
        
       x1=XB(k); x2=XB(k+1);
       y1=YB(k); y2=YB(k+1);

       u=((xo-x1)*(x2-x1)+(yo-y1)*(y2-y1))/((x2-x1)^2+(y2-y1)^2);

       if u<=0
           e(k)=(xo-x1)^2+(yo-y1)^2;
       elseif u>=1
           e(k)=(xo-x2)^2+(yo-y2)^2;
       else
           e(k)=(xo-x1-u*(x2-x1))^2+(yo-y1-u*(y2-y1))^2;
       end
       
    end

    x1=XB(LB); x2=XB(1);
    y1=YB(LB); y2=YB(1);
    
    u=((xo-x1)*(x2-x1)+(yo-y1)*(y2-y1))/((x2-x1)^2+(y2-y1)^2);

    if u<=0
        e(LB)=(xo-x1)^2+(yo-y1)^2;
    elseif u>=1
        e(LB)=(xo-x2)^2+(yo-y2)^2;
    else
        e(LB)=(xo-x1-u*(x2-x1))^2+(yo-y1-u*(y2-y1))^2;
    end
       
    db = sqrt(min(e));
    if db<1e-12
       db=0;
    end	

