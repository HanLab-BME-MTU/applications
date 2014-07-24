function distVector = calcContourDist(contourIn)

%
%Calculates arc length along a curve that may contain missing observations
%
%
%Hunter Elliott

%Determine number of points in contour
[chk,nPoints] = size(contourIn);

%Disable NaN warning
warning('off','MATLAB:chckxy:IgnoreNaN')

if chk ~=2
    disp('Input contour should be 2xM vector of X-Y coordinates!')
    return
end

if nPoints > 2
    
    if any(isnan(contourIn(:)))
    
        %Interpolate missing values
        contSpline = spline(1:nPoints,contourIn);
        intCurve = fnval(contSpline,1:nPoints);

    else
        intCurve = contourIn;
    end
    %Initialize distance vector
    distVector = zeros(nPoints,1);

    for j = 2:nPoints

        %Calculate distance from last point and add
        distVector(j) = distVector(j-1) + sqrt( (intCurve(1,j) - intCurve(1,j-1))^2+(intCurve(2,j) - intCurve(2,j-1))^2);


    end

elseif nPoints == 2
    %If only two points, just return the distance between them.
    distVector = norm(contourIn(:,1)-contourIn(:,2));
    
else
    distVector = 0;
end