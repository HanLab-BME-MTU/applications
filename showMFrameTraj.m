function showMFrameTraj(MFT,color,scale)
%showMFrameTraj Display the trajectories of speckles over multiframes.

figure(gcf);

base = MFT(:,1:2);
plot(base(:,2),base(:,1),[color '.']);
for j = 1:2:size(MFT,2)-3
   dispV = (MFT(:,j+2:j+3)-MFT(:,j:j+1))*scale;
   quiver(base(:,2), base(:,1), dispV(:,2), dispV(:,1), 0, color);
   base = base+dispV;
end

hold off;

