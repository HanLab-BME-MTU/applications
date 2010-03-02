function displaySpeckles(hAxes, tag, fileName, layerColor)

load(fileName);

if ~exist('cands', 'var')
    error('cands variable cannot be found in specified file.');
end

status = vertcat(cands(:).status); %#ok<NODEF>
cands = vertcat(cands(status == 1).Lmax);

line(cands(:,2),cands(:,1),'LineStyle', 'none', 'Marker', '.',...
    'Color',layerColor,'MarkerSize',6, 'Tag', tag, 'Parent', hAxes);
