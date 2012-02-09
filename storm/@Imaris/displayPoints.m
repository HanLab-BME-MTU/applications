function displayPoints(obj,points,spotSize,color,name)
if obj.displayEnabled
    % Create time data
    spotTime = zeros(size(points,1),1);
    
    % Set spot size
    if nargin > 2
        if numel(spotSize) == 1
            spotRadii = spotSize*ones(size(points,1),1);
        else
            spotRadii = spotSize;
        end
    else
        spotRadii = 5*ones(size(points,1),1);
    end
    
    % Create Imaris spot object
    imarisSpots = obj.imarisApp.mFactory.CreateSpots;
    imarisSpots.Set(points,spotTime,spotRadii);
    
    % Set spot color
    if nargin > 3
        imarisSpots.SetColor(color(1),color(2),color(3),color(4));
    else
        imarisSpots.SetColor(1.0,0.0,0.0,0.0);
    end
    
    % Name displayed in Imaris for this spot object
    if nargin > 4
        imarisSpots.mName = name;
    else
        imarisSpots.mName = 'Imaris: Points';
    end
    
    % Add the spot object to the scene
    obj.imarisApp.mSurpassScene.AddChild(imarisSpots);
end
end