function on_off_list=isOnCurve(pnts,curve)
% This function finds points on a pixel-position curve.

minx=min(vertcat(pnts(:,1),curve(:,1)));
miny=min(vertcat(pnts(:,2),curve(:,2)));

pnts(:,1)=pnts(:,1)-minx+1;
pnts(:,2)=pnts(:,2)-miny+1;
curve(:,1)=curve(:,1)-minx+1;
curve(:,2)=curve(:,2)-miny+1;

maxx=max(vertcat(pnts(:,1),curve(:,1)));
maxy=max(vertcat(pnts(:,2),curve(:,2)));

I=false(maxy,maxx);
idCurve = sub2ind([maxy,maxx],curve(:,2),curve(:,1));
idPnts  = sub2ind([maxy,maxx], pnts(:,2), pnts(:,1));

% Fill the matrix with the curve values:
I(idCurve)=true;

% Read out the values for the points of interest:
on_off_list=I(idPnts);