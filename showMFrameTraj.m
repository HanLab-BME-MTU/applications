function showMFrameTraj(MFT,index,color,scale)

figure(gcf);

base = MFT{index}(:,1:2);
plot(base(:,2),base(:,1),[color '.']);
for j = 1:2:size(MFT{index},2)-3
   dispV = (MFT{index}(:,j+2:j+3)-MFT{index}(:,j:j+1))*scale;
   quiver(base(:,2), base(:,1), dispV(:,2), dispV(:,1), 0, color);
   base = base+dispV;
end

hold off;

