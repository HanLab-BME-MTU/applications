function [p,t,dt,xNew,yNew,isSq]=delaunaySqLattice(xvec,yvec)
% first check if the given vectors are really representing a square
% lattice:
xPos=unique(xvec);
yPos=unique(yvec);

% size of the grid:
minX=min(xPos);
maxX=max(xPos);
minY=min(yPos);
maxY=max(yPos);

% the spacing between nodes:
dx=(maxX-minX)/(numel(xPos)-1);
dy=(maxY-minY)/(numel(yPos)-1);

% normalize xvec and yvec, to run from 1:numPts along x/y-dimension:
xvecNorm=(xvec-minX)/dx+1; % +1 for sub2ind
yvecNorm=(yvec-minY)/dy+1; % +1 for sub2ind

% the new normalized bounds 
maxXNorm=(maxX-minX)/dx+1;
maxYNorm=(maxY-minY)/dy+1;

% the position along the x/y-dimension:
xvecNew=linspace(1,maxXNorm,numel(xPos));
yvecNew=linspace(1,maxYNorm,numel(yPos));

% create the square lattice:
[X,Y] = meshgrid(xvecNew, yvecNew);
xNew=X(:);
yNew=Y(:);

% Check if all generated points agree with the input points by rewriting
% them in matrix form:
indOld=sub2ind([maxYNorm maxXNorm],yvecNorm,xvecNorm);
indNew=sub2ind([maxYNorm maxXNorm],yNew,xNew);

matOld=zeros(round(maxYNorm),round(maxXNorm));
matNew=zeros(round(maxYNorm),round(maxXNorm));

matOld(round(indOld))=1;
matNew(round(indNew))=1;

matDif=matOld-matNew;
isSq=~sum(matDif(:));

if ~isSq
    p=[];
    t=[];
    dt=[];
    xNew=xvec;
    yNew=yvec;
    return;
end

% Now we are sure that the input lattice is square. Generate the triangles
% now:

% the upper triangles:
pt1x =X(1:end-1,1:end-1); pt1x=pt1x(:);  
pt1y =Y(1:end-1,1:end-1); pt1y=pt1y(:);
pt2x =X(2:end  ,2:end);   pt2x=pt2x(:);
pt2y =Y(2:end  ,2:end);   pt2y=pt2y(:);
pt3x =X(1:end-1,2:end);   pt3x=pt3x(:);
pt3y =Y(1:end-1,2:end);   pt3y=pt3y(:);

idPt1=sub2ind([maxYNorm maxXNorm],pt1y,pt1x);
idPt2=sub2ind([maxYNorm maxXNorm],pt2y,pt2x);
idPt3=sub2ind([maxYNorm maxXNorm],pt3y,pt3x);

tupper=[idPt1 idPt2 idPt3];

% the lower triangles, pt1=pt1,pt3=pt2 as above:
idPt3=idPt2;
pt2x =X(2:end,1:end-1);   pt2x=pt2x(:);
pt2y =Y(2:end,1:end-1);   pt2y=pt2y(:);

idPt2=sub2ind([maxYNorm maxXNorm],pt2y,pt2x);

tlower=[idPt1 idPt2 idPt3];

% fuse all triangles:
t=vertcat(tupper,tlower);

%transform back:
xNew=(xNew-1)*dx+minX;
yNew=(yNew-1)*dy+minY;
p=[xNew yNew];

% triplot(t,xNew,yNew);
% hold on;
% plot(xvec,yvec,'or')
% hold off;

dt.X=p;
dt.Triangulation=t;

