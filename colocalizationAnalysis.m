% Francois Aguet, June 2010

function colocalizationAnalysis(rPath, gPath, fName, mode)

if nargin<4
    mode = 'WT';
end

chR = double(imread(rPath));
chG = double(imread(gPath));

if strcmp(mode, 'WT')
    % detection in both channels
    infoR = detectSpotsWT(chR);
    infoG = detectSpotsWT(chG);
    
    P1 = [infoR.xcom infoR.ycom];
    P2 = [infoG.xcom infoG.ycom];
elseif strcmp(mode, 'WTmerge')
    % merge channels
    % scale to [0..1]
    rangeR = max(chR(:))-min(chR(:));
    rangeG = max(chG(:))-min(chG(:));
    chRx = (chR-min(chR(:)))/rangeR;
    chGx = (chG-min(chG(:)))/rangeG;
    
    merge = max(max(chG(:)),max(chR(:)))*(chRx+chGx)/2; 
    infoM = detectSpotsWT(merge);
    P1 = [infoM.xcom infoM.ycom];
    P2 = P1;
else
    % red
    h = filterLoG(chR, 1.2);
    h(h<0) = 0;
    h = (h-min(h(:)))/(max(h(:))-min(h(:)));
    
    % local maxima
    h = locmax2d(h, [11 11]);
    maxh = h(h>0);
    
    L = graythresh(maxh);
    L = L/6;
    h(h<L) = 0;
        
    [yc xc] = find(h>0);
    P1 = [xc(:) yc(:)];
    
    % green
    h = filterLoG(chG, 1.2);
    h(h<0) = 0;
    h = (h-min(h(:)))/(max(h(:))-min(h(:)));
    
    % local maxima
    h = locmax2d(h, [11 11]);
    maxh = h(h>0);
    
    L = graythresh(maxh);
    L = L/6;
    h(h<L) = 0;
        
    [yc xc] = find(h>0);
    P2 = [xc(:) yc(:)];
end
    
    
m = size(P1,1);
n = size(P2,1);
dM = distMat2(P1,P2);
% search radius
nonLinkMarker = -1;
maxSR = 2;
m12 = nonLinkMarker*ones(m,m);
diagIdx = 1 + (0:m-1)*(m+1);
m12(diagIdx) = maxSR;
m21 = nonLinkMarker*ones(n,n);
diagIdx = 1 + (0:n-1)*(n+1);
m21(diagIdx) = maxSR;
m22 = ones(n,m);
M = [dM m12; m21 m22];

links12 = lap(M, -1);

map12 = [(1:m)' links12(1:m)];
map12(map12(:,2)>n,:) = [];

uniqueR = setdiff((1:m), map12(:,1));
uniqueG = setdiff((1:n), map12(:,2));


% display location of WT localized spots
figure;
imagesc(ch2rgb(chR, chG, [])); axis image;
%print('-depsc2', '-r300', [fName '_merge.eps']);
hold on;
plot(P1(:,1), P1(:,2), 'ro');
plot(P2(:,1), P2(:,2), 'gx');
plot([P1(map12(:,1),1)'; P2(map12(:,2),1)'], [P1(map12(:,1),2)'; P2(map12(:,2),2)'], 'b-');
%print('-depsc2', '-r300', [fName '_detectionAssignment.eps']);




P1x = P1(map12(:,1),1)';
P2x = P2(map12(:,2),1)';
P1y = P1(map12(:,1),2)';
P2y = P2(map12(:,2),2)';


sigma = 1.6;
w = ceil(2.5*sigma);
% loop through detected patches and fit Gaussians


% red channel
[ny nx] = size(chR);
xi = round(P1x);
yi = round(P1y);

% exclude detections at image border
idx = round(P1x)-w < 1 | round(P1x)+w > nx | round(P1y)-w < 1 | round(P1y)+w > ny |...
      round(P2x)-w < 1 | round(P2x)+w > nx | round(P2y)-w < 1 | round(P2y)+w > ny;

xi(idx) = [];
yi(idx) = [];

nc = length(xi);
xr = zeros(1,nc);
yr = zeros(1,nc);
Ar = zeros(1,nc);
cr = zeros(1,nc);
for k = 1:nc    
    window = chR(yi(k)-w:yi(k)+w, xi(k)-w:xi(k)+w);
    [p] = fitGaussian2D(window, [0 0 max(window(:)) sigma min(window(:))], 'xyAc');
    xr(k) = xi(k)+p(1);
    yr(k) = yi(k)+p(2);
    Ar(k) = p(3);
    cr(k) = p(5);
end


% green channel
xi = round(P2x);
yi = round(P2y);

xi(idx) = [];
yi(idx) = [];

nc = length(xi);
xg = zeros(1,nc);
yg = zeros(1,nc);
Ag = zeros(1,nc);
cg = zeros(1,nc);
for k = 1:nc    
    window = chG(yi(k)-w:yi(k)+w, xi(k)-w:xi(k)+w);
    [p] = fitGaussian2D(window, [0 0 max(window(:)) sigma min(window(:))], 'xyAc');
    xg(k) = xi(k)+p(1);
    yg(k) = yi(k)+p(2);
    Ag(k) = p(3);
    cg(k) = p(5);
end


% localize red channel only detections, read out in green
P1x = P1(uniqueR,1)';
P1y = P1(uniqueR,2)';
xi = round(P1x);
yi = round(P1y);
idx = round(P1x)-w < 1 | round(P1x)+w > nx | round(P1y)-w < 1 | round(P1y)+w > ny;
xi(idx) = [];
yi(idx) = [];

nc = length(xi);
xr_Rm = zeros(1,nc);
yr_Rm = zeros(1,nc);
Ar_Rm = zeros(1,nc);
Ag_Rm = zeros(1,nc);
% cg_Rm = zeros(1
for k = 1:nc    
    window = chR(yi(k)-w:yi(k)+w, xi(k)-w:xi(k)+w);
    [p] = fitGaussian2D(window, [0 0 max(window(:)) sigma min(window(:))], 'xyAc');
    xr_Rm(k) = xi(k)+p(1);
    yr_Rm(k) = yi(k)+p(2);
    Ar_Rm(k) = p(3);
    window = chG(yi(k)-w:yi(k)+w, xi(k)-w:xi(k)+w);
    [p] = fitGaussian2D(window, [p(1) p(2) max(window(:)) sigma min(window(:))], 'Ac');
    Ag_Rm(k) = p(3);
end



% localize green channel only detections, read out in red
P2x = P2(uniqueG,1)';
P2y = P2(uniqueG,2)';
xi = round(P2x);
yi = round(P2y);
idx = round(P2x)-w < 1 | round(P2x)+w > nx | round(P2y)-w < 1 | round(P2y)+w > ny;
xi(idx) = [];
yi(idx) = [];

nc = length(xi);
xg_Gm = zeros(1,nc);
yg_Gm = zeros(1,nc);
Ar_Gm = zeros(1,nc);
Ag_Gm = zeros(1,nc);
for k = 1:nc    
    window = chG(yi(k)-w:yi(k)+w, xi(k)-w:xi(k)+w);
    [p] = fitGaussian2D(window, [0 0 max(window(:)) sigma min(window(:))], 'xyAc');
    xg_Gm(k) = xi(k)+p(1);
    yg_Gm(k) = yi(k)+p(2);
    Ar_Gm(k) = p(3);
    window = chR(yi(k)-w:yi(k)+w, xi(k)-w:xi(k)+w);
    [p] = fitGaussian2D(window, [p(1) p(2) max(window(:)) sigma min(window(:))], 'Ac');
    Ag_Gm(k) = p(3);
end


% remove pairs with large distance after localization
d = sqrt((xr-xg).^2 + (yr-yg).^2);
idx = d>2;
xr(idx) = [];
yr(idx) = [];
xg(idx) = [];
yg(idx) = [];
Ar(idx) = [];
Ag(idx) = [];


% merge all detected pairs
xr = [xr xr_Rm xg_Gm];
yr = [yr yr_Rm yg_Gm];
Ar = [Ar Ar_Rm Ar_Gm];

xg = [xg xr_Rm xg_Gm];
yg = [yg yr_Rm yg_Gm];
Ag = [Ag Ag_Rm Ag_Gm];

fprintf('Red only: %d, green only: %d, total: %d\n', length(xr_Rm), length(xg_Gm), length(xr));



% idx = Ar > 500 | Ag > 500; % remove outliers
% xr(idx) = [];
% yr(idx) = [];
% xg(idx) = [];
% yg(idx) = [];
% Ar(idx) = [];
% Ag(idx) = [];



figure; imagesc(ch2rgb(chR, chG, [])); axis image;
hold on;
plot(xr, yr, 'ro');
plot(xg, yg, 'gx');
%print('-depsc2', '-r300', [fName '_localization.eps']);



ma = max(max(Ar), max(Ag));
R = colocalizationScatterPlot(Ar, Ag, [0 ma 0 ma], 'Red', 'Green');
title(['R^2 = ' num2str(R) ', C_{TfnR} = ' num2str(round(6*(1+sqrt(R))/(1-sqrt(R)))) ', N = ' num2str(length(xr))], 'FontName', 'Helvetica', 'FontSize', 16);% ', N = ' num2str(length(x))]);
if ~isempty(fName)
    print('-depsc2', '-r300', [fName '_scatterPlot.eps']);
end