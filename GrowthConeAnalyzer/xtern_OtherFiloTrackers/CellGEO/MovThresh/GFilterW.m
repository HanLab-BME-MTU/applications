function [X,Y]=GFilterW(x,y,n)
 
    if mod(n,2)==0
        n=n+1;
    end
    ots=floor(n/2);
    
    gs=(x(n)-x(1))/6;
    gx = linspace(-3*gs,3*gs,n);
    gf = exp( -(gx.^2)/(2*gs^2) );
    gf = gf / sum(gf);

    
    yi=(y(1)+y(2))/2;
    yf=(y(end-1)+y(end))/2;
    
    X=x; Y=X;
    yy=[yi*ones(1,ots) y yf*ones(1,ots)];
    
    for j=1:length(x)
        Y(j)=sum(yy(j:(j+n-1)).*gf);
    end

