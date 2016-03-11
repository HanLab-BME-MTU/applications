function [Th St M X Y]=AutoHist(F,Smoo)

    n=150;
    [hy,hx]=hist(F(:),n);

    [X,Y]=GFilterW(hx,hy,Smoo);
        
    yu=Y(3:end);
    yl=Y(1:(end-2));
    miind=find((Y(2:(end-1))-yu<0) & (Y(2:(end-1))-yl)<0);
    maind=find((Y(2:(end-1))-yu>0) & (Y(2:(end-1))-yl)>0);
    
    if isempty(miind)
        Th=0;
        St=0;
        M=0;
    else
              
        Xi=X(miind+1);
        Yi=Y(miind+1);
        Xa=X(maind+1);
        Ya=Y(maind+1);
        
        [maY,maI]=max(Ya);
        
        if Y(1)>maY           
            maX=X(1);
            maY=Y(1);
        else
            maX=Xa(maI);
        end
        
        miI=find(Xi>maX & Yi<maY/2);
        
        Th=Xi(miI(1));
        M=max(hy(hx>Th));
        St=1;
   
    end
    
   
    