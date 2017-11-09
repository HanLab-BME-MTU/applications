function nucleiPlot(nucleiStruc, nucNum, channelName)
%nucleiPlot Plots 3D half-transparent nucleus and takes 360 degreen
%snapshots
%   Detailed explanation goes here

close all
nucleus = nucleiStruc(nucNum).dapi;
v = nucleus;

% Make v a cube with zeros
if size(v,1) < size(v,2)
    v(size(v,1)+1:size(v,2),:,:)=0;
else if size(v,1) > size(v,2)
        v(:,size(v,2)+1:size(v,1),:)=0;
    end
end

[x,y,z]=meshgrid(1:size(v,1),1:size(v,2),1:size(v,3));

figure,
p = patch(isosurface(x,y,z,v));
% why not -3? see doc isosurface
isonormals(x,y,z,v,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1,1,39/152])
view(3); axis tight
camlight
lighting gouraud
% set white background
hold on
set(gcf,'color','white')
set(gca,'visible','off')
p.FaceAlpha = 0.382;

switch channelName
    case 'red'
        spots = nucleiStruc(nucNum).redSpot;
    case 'green'
        spots = nucleiStruc(nucNum).greenSpot;
    otherwise
        error('Channel not existed')
end

for i = 1:numel(spots)
    cordX = spots(i).cord(1);
    cordY = spots(i).cord(2);
    cordZ = spots(i).cord(3);
    plot3(cordX, cordY, cordZ, 'b', 'marker', 'o', 'markersize', 10, 'markerfacecolor', 'b');
end

light_handle = camlight('left');
rotSavePath = pwd;
% saveMeshRotation - rotate and save the current figure (could be used after plotMeshFigure)
 
% rotate the figure
for ct = 1:360
    camorbit(1,0,'camera')
    light_handle = camlight(light_handle,'headlight'); lighting phong
    drawnow
    
    % Save files
    toName = sprintf('rotate%03d',ct);
    saveas(gcf,fullfile(rotSavePath,toName), 'tiffn');
end

end

