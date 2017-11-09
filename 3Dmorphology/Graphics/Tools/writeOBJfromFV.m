function writeOBJfromFV(surface, savePath)


OBJ.vertices = surface.vertices;
OBJ.objects.data.vertices=surface.faces;
OBJ.objects.type='f';
write_wobj(OBJ, savePath);