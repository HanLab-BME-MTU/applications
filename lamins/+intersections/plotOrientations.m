function [ h ] = plotOrientations( g, numOrientations, scale, filter)
%plotOrientations Plot steerable filter orientations
%
% numOrientations - number of orientations to plot, 2 or 3?
% scale - scaling factor for the quiver vectors
% filter - binary filter mask of where to plot arrows... e.g. nms > 0

la = g.l.*(g.a-repmat(min(g.a,[],3),1,1,size(g.a,3)));
[M,MI] = sort(la,3);
N = size(M,3);

if(nargin < 2)
    numOrientations = N;
end
if(nargin < 3)
    scale = 5;
end
if(nargin < 4)
    filter = true([size(M,1) size(M,2)]);
end

if(numOrientations > 6)
    colors = parula(numOrientations);
else
    colors = [ 1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1];
end

[X,Y] = meshgrid(1:size(M,1),1:size(M,2));

hold on;
for i=0:numOrientations-1
    U = cos((MI(:,:,end-i)-1)*pi/N+pi/2).*M(:,:,end-i)*scale;
    V = sin((MI(:,:,end-i)-1)*pi/N+pi/2).*M(:,:,end-i)*scale;
%     U = cos(MI(:,:,end-i)*pi/N).*(i+1)*scale;
%     V = sin(MI(:,:,end-i)*pi/N).*(i+1)*scale;
    h(i+1) = quiver(X(filter),Y(filter),U(filter),V(filter),0,'color',colors(i+1,:),'AutoScale','off');
%     quiver(X(filter),Y(filter),U(filter),V(filter),0);
end
hold off;


end

