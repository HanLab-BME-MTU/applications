function saveMeshRotation(light_handle, rotSavePath)
 
% saveMeshRotation - rotate and save the current figure (could be used after plotMeshFigure)
 
% rotate the figure
for v=1:360
    camorbit(1,0,'camera')
    light_handle = camlight(light_handle,'headlight'); lighting phong
    drawnow
    toName = sprintf('rotate%03d',v);
    saveas(gcf,fullfile(rotSavePath,toName), 'tiffn');
end
