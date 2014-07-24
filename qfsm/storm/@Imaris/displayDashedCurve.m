function displayDashedCurve(obj,points,name)
if obj.displayEnabled
    % Number of points
    nPoints = size(points,1);
    
    % Create time data
    spotTime = zeros(nPoints,1);
    
    % Set spot size
    spotRadii = zeros(nPoints,1);
    
    % Create Imaris spot object
    imarisSpots = obj.imarisApp.mFactory.CreateSpots;
    imarisSpots2 = obj.imarisApp.mFactory.CreateSpots;
    imarisSpots.Set(points,spotTime,spotRadii);
    imarisSpots2.Set(points,spotTime,spotRadii);
    
    % Name displayed in Imaris for this spot object
    if nargin > 3
        imarisSpots.mName = name;
    else
        imarisSpots.mName = 'Imaris: Curve Points';
    end
    
    % Create edges
    edges1 = [(0:2:(nPoints-2))',(1:2:(nPoints-1))'];
    edges2 = [(1:2:(nPoints-1))',(2:2:nPoints)'];
    
    % Add the spot object to the scene
    imarisSpots.SetTrackEdges(edges1);
    imarisSpots2.SetTrackEdges(edges2);
    obj.imarisApp.mSurpassScene.AddChild(imarisSpots);
    obj.imarisApp.mSurpassScene.AddChild(imarisSpots2);
end
end