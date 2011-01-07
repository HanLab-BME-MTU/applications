function displayGrouping(hAxes, tag, G, layerColor)

V = G.V;
E = G.E;

i = E(:,1);
j = E(:,2);

X = [V(i,1) V(j,1)]';
Y = [V(i,2) V(j,2)]';

line(X(:), Y(:), 'Marker', 'none', 'Color', layerColor, 'Parent', hAxes, 'Tag', tag);
