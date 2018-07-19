maxZ = cell(4);
for i=1:4;
    maxZ{i} = max(cat(3,I{1,:}),[],3);
end;
maxZ_projection = [imadjust(maxZ{1}),imadjust(maxZ{2}),imadjust(maxZ{3}),imadjust(maxZ{4})];
imwrite(im2uint8(maxZ_projection),jet(256),'maxZprojection.png');
