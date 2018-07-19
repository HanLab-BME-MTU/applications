function [cx,cy,A]=myCentroid(x,y)

    rx=round(x);
    ry=round(y);
    n=max(rx)-min(rx)+1;
    m=max(ry)-min(ry)+1;
    sx=rx-min(rx)+1;
    sy=ry-min(ry)+1;
    BW=poly2mask(sx,sy,m,n);
    BW(m*(sx-1)+sy)=1;    
    [GX,GY]=meshgrid(1:n,1:m);
    MX=BW.*GX; MY=BW.*GY;
    A=sum(BW(:));
    cx=sum(MX(:))/A+min(rx)-1;
    cy=sum(MY(:))/A+min(ry)-1;
    
    
