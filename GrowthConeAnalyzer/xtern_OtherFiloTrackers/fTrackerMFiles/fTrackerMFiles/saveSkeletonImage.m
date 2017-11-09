function saveSkeletonImage(thisFrame, thisFrameSkeleton, thisSkeletonEnds, ...
    thisSkeletonJoints, thisCentroid, fileName, imagesFolder)

% Saves an image with the images and the skeletons overlayed

% Author: Santiago Costantino 
% santiago.costantino@umontreal.ca

radius=1;

green=adapthisteq(thisFrame);

blueEnds=zeros(size(thisFrame));
blueJoints=zeros(size(thisFrame));
for itEnds=1:size(thisSkeletonEnds, 1)
    thisBlueEnds=ImageRing(size(thisFrame), radius, thisSkeletonEnds(itEnds, 2),...
        thisSkeletonEnds(itEnds, 1), 1, 3);
    blueEnds=or(blueEnds, thisBlueEnds);
end
for itJoints=1:size(thisSkeletonJoints, 1)
    thisBlueJoints=ImageRing(size(thisFrame), radius, thisSkeletonJoints(itJoints, 2),...
        thisSkeletonJoints(itJoints, 1), 1, 3);
    blueJoints=or(blueJoints, thisBlueJoints);
end

centroidRing=ImageRing(size(thisFrame), 3*radius, round(thisCentroid(1,1)),round(thisCentroid(1,2)), 1, 5);

rgbFrame(:,:,2)=green;
rgbFrame(:,:,1)=or(thisFrameSkeleton, blueJoints)*255;
rgbFrame(:,:,3)=or(or(blueEnds, blueJoints),centroidRing)*255;

% imshow(rgbFrame);

imwrite(rgbFrame, [imagesFolder '\Skeletons\' fileName '.jpg'], 'jpg')
% hold on
% plot(thisSkeletonEnds(:,1),thisSkeletonEnds(:,2),'b*');
% plot(thisSkeletonJoints(:,1),thisSkeletonJoints(:,2),'y*');
% hold off


       
          