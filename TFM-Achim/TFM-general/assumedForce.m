function [force]=assumedForce(j,x,y)
xshift=1;
yshift=2;
if j==1
    force=x*0;
else
    force=(heaviside(1-((x-xshift).^2+(y-yshift).^2)).*(exp(-((x-xshift).^2+(y-yshift).^2))-1/exp(1)));
end

% if j==1
%     force=x*0;
% else
%     force=(heaviside(1-(x.^2+y.^2)).*(exp(-(x.^2+y.^2))-1/exp(1)));
% end