function dispOBJ(objPath)

OBJ = read_wobj(objPath);

FV.vertices = OBJ.vertices;
FV.faces = OBJ.objects.data.vertices;
colors = OBJ.vertices_texture(:,1);

figure
patch(FV, 'FaceColor', 'flat', 'cdata', colors); camlight
daspect([1 1 1]);