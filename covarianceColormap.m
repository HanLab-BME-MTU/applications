% Francois Aguet, October 2010

function covarianceColormap(M, paramNameList)

figure; imagesc(M); colormap(getMap()); axis image; caxis([-1 1]);
colorbar('YTick', -1:0.2:1);

np = size(M,1);
ta = 1:np;
set(gca, 'XTick', ta, 'YTick', ta, 'FontName', 'Helvetica', 'FontSize', 14);
if nargin>1 && length(paramNameList)==np
    set(gca, 'XTickLabel', paramNameList, 'YTickLabel', paramNameList);
else
    xlabel('Parameter index', 'FontName', 'Helvetica', 'FontSize', 14)
    ylabel('Parameter index', 'FontName', 'Helvetica', 'FontSize', 14)
end

M(eye(np)==1) = 0;
[x,y] = ind2sub(size(M), find(abs(M)>=0.8));
hold on;
plot(x,y, 'xy', 'LineWidth', 3, 'MarkerSize', 20);



function map = getMap()
values = -1:1/100:1;
N = length(values);
map = zeros(N,3);

ridx = values<0;
map(ridx,1) = -values(ridx);
gidx = values>0;
map(gidx,2) = values(gidx);