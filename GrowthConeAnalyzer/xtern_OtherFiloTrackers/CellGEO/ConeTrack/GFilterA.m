function Y=GFilterA(y,n)
 
    if length(y)<n || n<=1
        Y=y;
        return;
    end

    if mod(n,2)==0
        n=n+1;
    end
    ots=floor(n/2);
    
    gs=(n-1)/6;
    gx = linspace(-3*gs,3*gs,n);
    gf = exp( -(gx.^2)/(2*gs^2) );
    gf = gf / sum(gf);
    
    yi=(y(1)+y(2))/2;
    yf=(y(end-1)+y(end))/2;
    
    Y=y; 
    if size(y,1)==1
        yy=[yi*ones(1,ots) y yf*ones(1,ots)];
    elseif size(y,2)==1
        yy=[yi*ones(1,ots) y' yf*ones(1,ots)];
    else
        Y=[];
        return;
    end
            
    for j=1:length(y)
        Y(j)=sum(yy(j:(j+n-1)).*gf);
    end

