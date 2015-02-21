function plot_color_dot_from_matrix(inMatrix,inMarkersize, min_value,max_value)
% function to plot a matrix as a dot map with corresponding color, default
% as 'jet'

% Liya 2015.02

colorArray =colormap('jet');

inMatrix = double(inMatrix);

if(nargin<2)
    inMarkersize = 24;
end
if(nargin<3)
    min_value = min(min(inMatrix(~isnan(inMatrix))));
end
if(nargin<4)
    max_value = (max(max(inMatrix(~isnan(inMatrix)))) - median(inMatrix(~isnan(inMatrix))))*1.2 + median(inMatrix(~isnan(inMatrix)));
end

% normalize matrix to 0~1
inMatrix =  double(inMatrix-min_value)/(max_value-min_value);

% give the color index
colorIndexMatrix = round(inMatrix*size(colorArray,1));
colorIndexMatrix(colorIndexMatrix<1)=1;
colorIndexMatrix(colorIndexMatrix>size(colorArray,1))=size(colorArray,1);

% build meshgrid
[X,Y] = meshgrid(1:size(inMatrix,2),1:size(inMatrix,1));

%display the plot
hold off;
% show to set figure size
imagesc(inMatrix+nan);axis image;

curunits = get(gca, 'Units');
set(gca, 'Units', 'Points');
cursize = get(gca, 'Position');
% set the verticle to the same as horizontal
% if(size(inMatrix,1)>size(inMatrix,2))
cursize(3) = size(inMatrix,2)*inMarkersize;
cursize(4) = size(inMatrix,1)*inMarkersize;
set(gca, 'Position',cursize);
% end
cursize = get(gca, 'Position');

set(gcf, 'Units', 'Points');
set(gcf, 'Position',cursize+[0 0 72 72]);

set(gca, 'Units', curunits);
set(gcf, 'Units', curunits);

if(size(inMatrix,2)<size(inMatrix,1)/1.5)
    inMarkersize = inMarkersize/1.5;
end
% plot dot by dot
for i = 1 : numel(X)
    % dot
    if ~isnan(Y(i)) && ~isnan(colorIndexMatrix(i))
    plot(X(i)-0.5,Y(i)-0.5,'.','markersize',63*inMarkersize/24,'color',colorArray(colorIndexMatrix(i),:));
    hold on;
    % circles
    plot(X(i)-0.5,Y(i)-0.5,'bo','markersize',21*inMarkersize/24);
    end
end
axis equal;
axis([0 size(inMatrix,2)+2  0 size(inMatrix,1)+2]);
box off;
plot([0 size(inMatrix,2)],[size(inMatrix,1) size(inMatrix,1)]);
plot([size(inMatrix,2) size(inMatrix,2)],[0 size(inMatrix,1)]);
