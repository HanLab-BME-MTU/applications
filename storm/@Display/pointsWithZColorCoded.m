function pointsWithZColorCoded(obj)

colorTable = hsv(obj.nColorSlices);
zStepPerColor = (max(obj.data.points(:,3))-min(obj.data.points(:,3)))/obj.nColorSlices;

for color=1:obj.nColorSlices
    % Assign the points to a color
    pointIdxsSameColor = find(obj.data.points(:,3) >= min(obj.data.points(:,3))+(color-1)*zStepPerColor & obj.data.points(:,3) <= min(obj.data.points(:,3))+color*zStepPerColor);
    
    % Create Imaris Spot object
    if size(pointIdxsSameColor) > 0
        colorRGBA = [colorTable(color,1) colorTable(color,2) colorTable(color,3) 0.0];
        obj.imaris.displayPoints(obj.data.points(pointIdxsSameColor,:),3,colorRGBA,'Display: Points');
    end
end

obj.imaris.fitCamera();

end

