function rgbFrame = showSkeletonImage(thisFrame, thisFrameSkeleton, thisSkeletonEnds, thisSkeletonJoints,show)

% Saves an image with the images and the skeletons overlayed
% Author: Santiago Costantino 
% santiago.costantino@umontreal.ca


radius=1;

green=(thisFrame);

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

rgbFrame(:,:,2)=green;
rgbFrame(:,:,1)=or(thisFrameSkeleton, blueJoints)*255;
rgbFrame(:,:,3)=or(blueEnds, blueJoints)*255;

if show == 1
    figure('Name', 'Cone and Skeleton');imshow(rgbFrame);
end

% hold on
% plot(thisSkeletonEnds(:,1),thisSkeletonEnds(:,2),'b*');
% plot(thisSkeletonJoints(:,1),thisSkeletonJoints(:,2),'y*');
% hold off


       
          