% Francois Aguet, October 2010

function plotCorrelationMatrix(M, paramNameList)

triuMask = triu(ones(size(M)));
[x,y] = ind2sub(size(M), find(abs(M)>0.8 & triuMask==0));

% remove upper triangular part
M = M-triu(M);

% convert values [-1..1] to RGB
M = cat(3, -M.*(M<0), M.*(M>0), zeros(size(M)));

M(repmat(triuMask, [1 1 3])==1)=1;


figure('Position', [440 378 300 250]);
imagesc(M); colormap(getMap()); axis image; caxis([-1 1]);
colorbar('YTick', -1:0.2:1);

np = size(M,1);
ta = 1:np;
set(gca, 'XTick', ta, 'YTick', ta, 'XAxisLocation', 'top', 'TickLength', [0 0],...
    'FontName', 'Helvetica', 'FontSize', 14);
if nargin>1 && length(paramNameList)==np
    set(gca, 'XTickLabel', paramNameList, 'YTickLabel', paramNameList);
else
    xlabel('Parameter index', 'FontName', 'Helvetica', 'FontSize', 14)
    ylabel('Parameter index', 'FontName', 'Helvetica', 'FontSize', 14)
end
box off;

hold on;
plot(y, x, '.y', 'LineWidth', 3, 'MarkerSize', 20);



function map = getMap()
values = -1:1/100:1;
N = length(values);
map = zeros(N,3);

ridx = values<0;
map(ridx,1) = -values(ridx);
gidx = values>0;
map(gidx,2) = values(gidx);
