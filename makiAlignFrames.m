function dataStruct = makiAlignFrames(dataStruct)
%MAKIROTATEFRAMES aligns frames based on translation and rotation estimates derived from tracks
%
%SYNOPSIS dataStruct = makiAlignFrames(dataStruct)
%
%INPUT  dataStruct : dataStruct as in makiMakeDataStruct with at least the
%                    fields "dataProperties", "initCoord", "planeFit", 
%                    "tracks" and "updatedClass". 
%                    Optional. Loaded interactively if not input.
%
%OUTPUT dataStruct : Same as input, with added field "frameAlignment".
%                    For each frame, it contains the subfields:
%           .centerOfMass: Center of mass of features in frame.
%           .eulerAnglesX: The Euler angles defining the rotation of this
%                          frame from its reference frame (which is the 4th
%                          entry in eulerAnglesX). Euler angles follow the
%                          x-convention. NaN means that no rotation was
%                          estimated for this frame, its coordinate system
%                          is obtained via the plane fit.
%           .coordSystem : Coordinate system of frame. If eulerAnglesX are
%                          NaN, it comes from the plane fit. If
%                          eulerAnglesX are not NaN, it is obtained by
%                          rotating the coordinate system of the reference
%                          frame.
%           .alignedCoord: Aligned coordinates in each frame, where both
%                          center of mass translation and rigid body
%                          rotation are compensated for.
%
%Khuloud Jaqaman, July 2007

%% preamble

%load dataStruct if not input
if nargin == 0 || isempty(dataStruct)
    dataStruct = makiLoadDataFile;
end

%get number of frames in movie
numFrames = dataStruct.dataProperties.movieSize(end);

%put tracks into matrix format
[tracksInfo,tracksIndx] = convStruct2MatNoMS(dataStruct.tracks);

%get total number of tracks
numTracks = size(tracksIndx,1);

%assign 0 to outlier (unaligned & lagging) features in each frame
for iFrame = 1 : numFrames
    outIndx = [dataStruct.updatedClass(iFrame).unalignedIdx dataStruct.updatedClass(iFrame).laggingIdx];
    for iFeat = outIndx
        tracksIndx(tracksIndx(:,iFrame)==iFeat,iFrame) = 0;
    end
end

%reserve memory for euler angles and rotation frame of reference
eulerAnglesX = NaN(numFrames,4);

%find which frames need their rotation to be estimated
frames2align = [];
for t = 1 : numFrames
    if isempty(dataStruct.planeFit(t).planeVectors)
        frames2align = [frames2align t];
    end
end
framesOK = setxor((1:numFrames),frames2align);

%if all frames are empty, make the first frame non-empty (it will be taken
%as the first reference frame)
if isempty(framesOK)
    frames2align = frames2align(2:end);
    framesOK = 1;
end

%find empty frames that are before the first non-empty frame
framesBefore1 = frames2align(frames2align < framesOK(1));
framesAfter1 = setxor(frames2align,framesBefore1); % the rest of the frames

%the reference frame of each empty frame before the first non-empty one is
%the frame after it
eulerAnglesX(framesBefore1,4) = framesBefore1 + 1;

%for the rest, their reference frame is the frame before each of them
eulerAnglesX(framesAfter1,4) = framesAfter1 - 1;

%% center of mass shift

%get the center of mass of each frame
centerOfMass = vertcat(dataStruct.planeFit.planeOrigin);

%shift the track coordinates so that the origin is at the center of mass
%in each frame
for iFrame = 1 : numFrames
    tracksInfo(:,8*(iFrame-1)+1:8*(iFrame-1)+3) = ...
        tracksInfo(:,8*(iFrame-1)+1:8*(iFrame-1)+3) - ...
        repmat(centerOfMass(iFrame,:),numTracks,1);
end

%% estimation of rotation 

%go over all empty frames to estimate their rotation with respect to their
%reference frame
if ~isempty(frames2align)
    for iFrame = frames2align

        %determine reference frame
        jFrame = eulerAnglesX(iFrame,4);

        %get the indices of features in this frame
        featIndx1 = tracksIndx(:,iFrame);

        %get the indices of features in reference frame
        featIndx2 = tracksIndx(:,jFrame);

        %find those features that are linked between the two frames
        goodIndx = find(featIndx1 ~= 0 & featIndx2 ~= 0);

        %get their coordinates
        coord1 = tracksInfo(goodIndx,(jFrame-1)*8+1:(jFrame-1)*8+3);
        coord2 = tracksInfo(goodIndx,(iFrame-1)*8+1:(iFrame-1)*8+3);

        %estimate rotation angles (Euler angles following the x-convention)
        x0 = [0 0 0]; %initial guess
        lb = [0 0 0]; %lower bound
        ub = pi*[2 1 2]; %upper bound
        options = optimset('Jacobian','on','Display','off'); %minimization option - use analytical Jacobian
        rotationAngles = lsqnonlin(@calcRotateResiduals,x0,lb,ub,options,coord1,coord2);

        %save rotation angles for output
        eulerAnglesX(iFrame,1:3) = rotationAngles;

    end
end

%% rotation of frame coordinate systems

%reserve memory
coordSystem = repmat(eye(3),[1 1 numFrames]);

%get the coordinate systems of non-empty frames from the plane fit
for iFrame = framesOK
    if isempty(dataStruct.planeFit(iFrame).planeVectors)
        coordSystem(:,:,iFrame) = eye(3);
    else
        coordSystem(:,:,iFrame) = dataStruct.planeFit(iFrame).planeVectors;
    end
end    

%calculate the coordinate systems of empty frames that are before the first
%non-empty frame
%go from last to first
if ~isempty(framesBefore1)
    for iFrame = framesBefore1(end:-1:1)

        %get the reference frame
        jFrame = eulerAnglesX(iFrame,4);

        %fetch the Euler angles for this frame
        phi = eulerAnglesX(iFrame,1);
        theta = eulerAnglesX(iFrame,2);
        psi = eulerAnglesX(iFrame,3);

        %calculate the rotation matrix
        rotationMat1 = [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1];
        rotationMat2 = [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];
        rotationMat3 = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];
        rotationMat = rotationMat1 * rotationMat2 * rotationMat3;

        %rotate the coordinate system of this frame
        coordSystem(:,:,iFrame) = rotationMat * coordSystem(:,:,jFrame);

    end
end

%calculate the coordinate systems of empty frames that are after the first
%non-empty frame
%go from first to last
if ~isempty(framesAfter1)
    for iFrame = framesAfter1

        %get the reference frame
        jFrame = eulerAnglesX(iFrame,4);

        %fetch the Euler angles for this frame
        phi = eulerAnglesX(iFrame,1);
        theta = eulerAnglesX(iFrame,2);
        psi = eulerAnglesX(iFrame,3);

        %calculate the rotation matrix
        rotationMat1 = [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1];
        rotationMat2 = [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];
        rotationMat3 = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];
        rotationMat = rotationMat1 * rotationMat2 * rotationMat3;

        %rotate the coordinate system of this frame
        coordSystem(:,:,iFrame) = rotationMat * coordSystem(:,:,jFrame);

    end
end

%% coordinates in rotated & translated coordinate system

%reserve memory
alignedCoord(1:numFrames,1) = struct('values',[]);

%go over all frames
for iFrame = 1 : numFrames

    %get feature coordinates in this frame
    coordTmp = dataStruct.initCoord(iFrame).allCoord;
    numFeatures = size(coordTmp,1);
    
    %subtract center of mass from coordinates
    coordTmp(:,1:3) =  coordTmp(:,1:3) - repmat(centerOfMass(iFrame,:),numFeatures,1);
    
    %get rotation matrix
    rotationMat = inv(coordSystem(:,:,iFrame));
    
    %calculate matrix to propagate error from original to rotated
    %coordinates
    errorPropMat = rotationMat.^2;
    
    %calculate coordinates in rotated coordinate system
    alignedCoord(iFrame).values(:,1:3) = (rotationMat * coordTmp(:,1:3)')';
    
    %propagate error
    alignedCoord(iFrame).values(:,4:6) = sqrt((errorPropMat*((coordTmp(:,4:6)).^2)')');

end

%% output to dataStruct

%fill in field
frameAlignment(1:numFrames,1) = struct('centerOfMass',[],'eulerAnglesX',[],...
    'coordSystem',[],'alignedCoord',[]);
for iFrame = 1 : numFrames
    frameAlignment(iFrame).centerOfMass = centerOfMass(iFrame,:);
    frameAlignment(iFrame).eulerAnglesX = eulerAnglesX(iFrame,:);
    frameAlignment(iFrame).coordSystem = coordSystem(:,:,iFrame);
    frameAlignment(iFrame).alignedCoord = alignedCoord(iFrame).values;
end

%write field into dataStruct
dataStruct.frameAlignment = frameAlignment;


%% subfunction that calculates rotation residuals
function [F,J] = calcRotateResiduals(x0,coord1,coord2)

%fetch the Euler angles
phi = x0(1);
theta = x0(2);
psi = x0(3);

%calculate the rotation matrix
rotationMatPsi = [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1];
rotationMatTheta = [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];
rotationMatPhi = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];
rotationMat = rotationMatPsi * rotationMatTheta * rotationMatPhi;

%calculate the residuals and put them in a 1D vector
F = coord2 - coord1 * rotationMat';
F = F(:);

%calculate the Jacobian matrix - initialize
J = zeros(length(F),3);

%calculate the derivative of the rotation matrix with respect to phi
rotationMatPhiD = [-sin(phi) cos(phi) 0; -cos(phi) -sin(phi) 0; 0 0 0];
rotationMat = rotationMatPsi * rotationMatTheta * rotationMatPhiD;

%calculate first column of Jacobian matrix
JacobianCol = -coord1 * rotationMat';
J(:,1) = JacobianCol(:);

%calculate the derivative of the rotation matrix with respect to theta
rotationMatThetaD = [0 0 0; 0 -sin(theta) cos(theta); 0 -cos(theta) -sin(theta)];
rotationMat = rotationMatPsi * rotationMatThetaD * rotationMatPhi;

%calculate second column of Jacobian matrix
JacobianCol = -coord1 * rotationMat';
J(:,2) = JacobianCol(:);

%calculate the derivative of the rotation matrix with respect to psi
rotationMatPsiD = [-sin(psi) cos(psi) 0; -cos(psi) -sin(psi) 0; 0 0 0];
rotationMat = rotationMatPsiD * rotationMatTheta * rotationMatPhi;

%calculate third column of Jacobian matrix
JacobianCol = -coord1 * rotationMat';
J(:,3) = JacobianCol(:);

