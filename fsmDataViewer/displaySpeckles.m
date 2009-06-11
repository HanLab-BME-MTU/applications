function displaySpeckles(hAxes, tag, fileName, layerColor)

load(fileName);

if ~exist('cands', 'var')
    error('cands variable cannot be found in specified file.');
end

% Extract speckle classes
primary = find([cands.status] == 1 & [cands.speckleType] == 1);
secondary = find([cands.status] == 1 & [cands.speckleType] == 2);
tertiary = find([cands.status] == 1 & [cands.speckleType] == 3);
higher = find([cands.status] == 1 & [cands.speckleType] > 3);

% Extract speckle positions
pPos=reshape([cands(primary).Lmax], 2, length([cands(primary).Lmax])/2)';
sPos=reshape([cands(secondary).Lmax], 2, length([cands(secondary).Lmax])/2)';
tPos=reshape([cands(tertiary).Lmax], 2, length([cands(tertiary).Lmax])/2)';
hPos=reshape([cands(higher).Lmax], 2, length([cands(higher).Lmax])/2)';

% Primary speckles
line(pPos(:,2),pPos(:,1),'LineStyle', 'none',...
    'Marker', '.','Color',layerColor,'MarkerSize',6,...
    'Tag', tag, 'Parent', hAxes);

% Secondary peckles
line(sPos(:,2),sPos(:,1),'LineStyle', 'none',...
    'Marker', '+','Color',layerColor,'MarkerSize',6,...
    'Tag', tag, 'Parent', hAxes);

% Tertiary speckles
line(tPos(:,2),tPos(:,1),'LineStyle', 'none',...
    'Marker', '^','Color',layerColor,'MarkerSize',6,...
    'Tag', tag, 'Parent', hAxes);

% Higher-order speckles
line(hPos(:,2),hPos(:,1),'LineStyle', 'none',...
    'Marker', '*','Color',layerColor,'MarkerSize',6,...
    'Tag', tag, 'Parent', hAxes);

end