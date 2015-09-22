function displayVectors(obj,vectors)
% The vectors are displayed at the origin
pointsStart = zeros(size(vectors));
pointsEnd = vectors;
lenVectors = sqrt(sum(vectors.^2,2));

% Display the vectors
nameSegments = 'Imaris: Vectors';
obj.displaySegments(pointsStart,pointsEnd,nameSegments);

% Display the vector tips
nameVectorTips = 'Imaris: Vectors Tips';
obj.displayPoints(pointsEnd,lenVectors/10,[1.0 0.0 0.0 0.0],nameVectorTips);
end