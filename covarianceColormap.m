% Francois Aguet, October 2010

function covarianceColormap(M)

figure; imagesc(M); colormap(getMap()); axis image; colorbar; caxis([-1 1]);
xlabel('Rate constant index', 'FontName', 'Helvetica', 'FontSize', 14)
ylabel('Rate constant index', 'FontName', 'Helvetica', 'FontSize', 14)



function map = getMap()
% values = sort(values);
values = -1:1/100:1;
N = length(values);
map = zeros(N,3);

ridx = values<0;
map(ridx,1) = -values(ridx);
gidx = values>0;
map(gidx,2) = values(gidx);