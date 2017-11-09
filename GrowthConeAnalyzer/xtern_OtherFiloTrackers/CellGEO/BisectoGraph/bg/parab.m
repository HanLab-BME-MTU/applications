function [X,Y]=parab(t,xf,yf,xr,yr,ph)

    x=xr+abs(t)*cos(ph);
    y=yr+abs(t)*sin(ph);

    x0=x; y0=y; an=ph+pi/2;

    A1=-sin(an); B1= cos(an);
    C1= x0*sin(an)-y0*cos(an);

    tm=atan2(y-yf,x-xf);

    x0=(x+xf)/2; y0=(y+yf)/2; an=tm+pi/2;

    A2=-sin(an); B2= cos(an);
    C2= x0.*sin(an)-y0.*cos(an);

    X=(B1.*C2-B2.*C1)./(A1.*B2-A2.*B1);
    Y=(C1.*A2-C2.*A1)./(A1.*B2-A2.*B1);

