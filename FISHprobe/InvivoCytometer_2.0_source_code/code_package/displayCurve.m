function displayCurve(obj,points,name)
if obj.displayEnabled
    % Number of points
    nPoints = size(points,1);
    
    % Create time data
    spotTime = zeros(nPoints,1);
    
    % Set spot size
    spotRadii = zeros(nPoints,1);
    
    % Create Imaris spot object
    imarisSpots = obj.imarisApp.mFactory.CreateSpots;
    imarisSpots.Set(points,spotTime,spotRadii);
    
    % Name displayed in Imaris for this spot object
    if nargin > 3
        imarisSpots.mName = name;
    else
        imarisSpots.mName = 'Imaris: Curve Points';
    end
    
    % Add the spot object to the scene
    imarisSpots.SetTrackEdges([(1:(nPoints-1))',(2:nPoints)']-1);
    obj.imarisApp.mSurpassScene.AddChild(imarisSpots);
end
end