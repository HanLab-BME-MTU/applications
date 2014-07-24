function displaySegments(obj,pointsStart,pointsEnd,name)
if obj.displayEnabled
    % Number of segments
    nSeg = size(pointsStart,1);
    
    % Set spot size
    spotRadii = zeros(2*nSeg,1);
    
    % Create edges
    edges = [(0:nSeg-1)',(nSeg:(2*nSeg)-1)'];
    
    % Create points
    points = [pointsStart;pointsEnd];
    
    % Create time data
    spotTime = zeros(2*nSeg,1);
    
    % Create Imaris spot object
    imarisSpots = obj.imarisApp.mFactory.CreateSpots;
    imarisSpots.Set(points,spotTime,spotRadii);
    
    % Name displayed in Imaris for this spot object
    if nargin > 3
        imarisSpots.mName = name;
    else
        imarisSpots.mName = 'Imaris: Segments';
    end
    
    % Add the spot object to the scene
    imarisSpots.SetTrackEdges(edges);
    
    obj.imarisApp.mSurpassScene.AddChild(imarisSpots);
end
end