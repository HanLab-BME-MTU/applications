function bandMapPlot(image, bandMap, handles)
% close all
axes(handles.axes2);
% set(gcf,'position',get(0,'screensize'))
imshow(imcomplement(image),[]),
title('Detected Bands'), hold on
% p goes row by row and q goes column by column
for p = 1:size(image,1)
    for q = 1:size(image,2)
        if bandMap(p,q) == 1
            plot(q, p, 'r.'),hold on
        end
    end
end
hold off