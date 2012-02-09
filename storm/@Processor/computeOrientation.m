function computeOrientation(obj,filterLength,filterAngularSampling)

meanLocPrec = mean(obj.data.error,1);

obj_data_orientation = zeros(obj.data.nPoints,3);
obj_data_magnitude = zeros(obj.data.nPoints,1);

radii = repmat(filterLength/2,obj.data.nPoints,1);

[idxs,~] = KDTreeBallQuery(obj.data.points,obj.data.points,radii);

phiAngles = 0:filterAngularSampling:180-filterAngularSampling;
thetaAngles = 0:filterAngularSampling:180-filterAngularSampling;
[phi,theta] = meshgrid(phiAngles,thetaAngles);

obj_data_nPoints = obj.data.nPoints;
obj_data_points = obj.data.points;

directionVector = zeros(1,3,size(theta,1),size(theta,2));
directionVector(1,1,:,:) = sind(theta).*cosd(phi)/meanLocPrec(1);
directionVector(1,2,:,:) = sind(theta).*sind(phi)/meanLocPrec(2);
directionVector(1,3,:,:) = cosd(theta)/meanLocPrec(3);

parfor i=1:obj_data_nPoints
    centerPoint = obj_data_points(i,:);
    centerPoint = centerPoint./meanLocPrec;
    
    directionVectorGrid = repmat(directionVector,[size(idxs{i},1),1,1,1]);
    
    points = obj_data_points(idxs{i},:);
    points = points./repmat(meanLocPrec,size(points,1),1);
    points = repmat(points,[1,1,size(phiAngles,2),size(thetaAngles,2)]);
    
    linePoints = repmat(centerPoint,[size(idxs{i},1),1,size(phiAngles,2),size(thetaAngles,2)]);
    
    % Compute distances
    t = sum(directionVectorGrid.*(points-linePoints),2)./sum(directionVectorGrid.*directionVectorGrid,2);
    
    tNorm = t/filterLength*2;
    classification = zeros(size(t));
    
    %     classification(tNorm>=0) = 1;
    %     classification(tNorm<0) = 2;
    
    %     classification(tNorm>=0.5) = 1;
    %     classification(tNorm>=0&tNorm<0.5) = 2;
    %     classification(tNorm<0&tNorm>-0.5) = 3;
    %     classification(tNorm<=-0.5) = 4;
    
    classification(tNorm>=0.75) = 1;
    classification(tNorm>=0.5&tNorm<0.75) = 2;
    classification(tNorm>=0.25&tNorm<0.5) = 3;
    classification(tNorm>=0&tNorm<0.25) = 4;
    classification(tNorm<0&tNorm>-0.25) = 5;
    classification(tNorm<=-0.25&tNorm>-0.5) = 6;
    classification(tNorm<=-0.5&tNorm>-0.75) = 7;
    classification(tNorm<=-0.75) = 8;
    
    t = repmat(t,[1,3,1,1]);
    orthogonalVectors = points-(linePoints+t.*directionVectorGrid);
    distanceSq = sum(orthogonalVectors.*orthogonalVectors,2); % Squared distances
    
    % Compute the exponential value of this distance
    contribution = exp(-distanceSq/2);
    
    minBinResponse = 10;
    
    %     response = (minBinResponse+sum(contribution.*(classification==1),1)) ...
    %         .*(minBinResponse+sum(contribution.*(classification==2),1));
    %
    %     response = (minBinResponse+sum(contribution.*(classification==1),1)) ...
    %         .*(minBinResponse+sum(contribution.*(classification==2),1)) ...
    %         .*(minBinResponse+sum(contribution.*(classification==3),1)) ...
    %         .*(minBinResponse+sum(contribution.*(classification==4),1));
    
    response = (minBinResponse+sum(contribution.*(classification==1),1)) ...
        .*(minBinResponse+sum(contribution.*(classification==2),1)) ...
        .*(minBinResponse+sum(contribution.*(classification==3),1)) ...
        .*(minBinResponse+sum(contribution.*(classification==4),1)) ...
        .*(minBinResponse+sum(contribution.*(classification==5),1)) ...
        .*(minBinResponse+sum(contribution.*(classification==6),1)) ...
        .*(minBinResponse+sum(contribution.*(classification==7),1)) ...
        .*(minBinResponse+sum(contribution.*(classification==8),1));
    
    %     response = (minBinResponse+sum(contribution.*(classification==1),1)) ...
    %         +(minBinResponse+sum(contribution.*(classification==2),1)) ...
    %         +(minBinResponse+sum(contribution.*(classification==3),1)) ...
    %         +(minBinResponse+sum(contribution.*(classification==4),1)) ...
    %         +(minBinResponse+sum(contribution.*(classification==5),1)) ...
    %         +(minBinResponse+sum(contribution.*(classification==6),1)) ...
    %         +(minBinResponse+sum(contribution.*(classification==7),1)) ...
    %         +(minBinResponse+sum(contribution.*(classification==8),1));
    
    response = squeeze(response);
    
    % Find the max response of all the different orientations
    [rowVal rowIdx] = max(response,[],1);
    [colVal colIdx] = max(rowVal);
    rowIdx = rowIdx(colIdx);
    phiMax = phiAngles(colIdx);
    thetaMax = thetaAngles(rowIdx);
    
    obj_data_orientation(i,:) = [sind(thetaMax)*cosd(phiMax) sind(thetaMax)*sind(phiMax) cosd(thetaMax)];
    obj_data_magnitude(i) = colVal;
end

obj.data.orientation = obj_data_orientation;
obj.data.magnitude = obj_data_magnitude;

disp('Process: Orientation computed!');

end

