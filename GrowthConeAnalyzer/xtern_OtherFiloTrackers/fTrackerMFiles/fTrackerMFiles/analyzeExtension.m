% Track the ends that are the result of analyzeSkeletonMovie and meassure
% the distance to the closest joint

% Author: Santiago Costantino 
% santiago.costantino@umontreal.ca


clear

imagesFolder='Images\Chris\Series6\';

%% Read parameters form settings file

settingsRead=inifile([imagesFolder 'Settings.ini'],'readall');
fileName=getFromINI(settingsRead, 'fileName');
maxDisplacement=getFromINI(settingsRead, 'maxDisplacement');
frames=getFromINI(settingsRead, 'frames');
trackingParameters.mem=getFromINI(settingsRead, 'lostFrames');
trackingParameters.dim=getFromINI(settingsRead, 'dimensions');
trackingParameters.good=getFromINI(settingsRead, 'minLenght');
trackingParameters.quiet=0;


fullEndInfo=[];

%% Load data
skeletonEnds=load([imagesFolder 'Tracks\' fileName 'Ends.dat']);
skeletonJoints=load([imagesFolder 'Tracks\' fileName 'Joints.dat']);
skeletonCentroids=load([imagesFolder 'Tracks\' fileName 'CenterOfMass.dat']);

%% Track
trackResult = track(skeletonEnds, maxDisplacement, trackingParameters);

%% Analize tracks

for itTracks=1:max(trackResult(:,4))
    for itEnds=find(trackResult(:,4)==itTracks)'
        thisEnd=[trackResult(itEnds,1), trackResult(itEnds,2)];
        theseJoints=[skeletonJoints(find(skeletonJoints(:,3)==trackResult(itEnds,3)),1),...
            skeletonJoints(find(skeletonJoints(:,3)==trackResult(itEnds,3)),2)];
        thisTime=max(skeletonJoints(find(skeletonJoints(:,3)==2),3));
        [thisDistance, closestJoint]=getDistance(thisEnd, theseJoints);
        [thisDistanceToCentoid, centroid]=getDistance(thisEnd, ...
            [skeletonCentroids(thisTime, 2), skeletonCentroids(thisTime, 1)]);
        fullEndInfo=[fullEndInfo; thisEnd closestJoint thisDistance...
            skeletonCentroids(trackResult(itEnds,3), 2) ...
            skeletonCentroids(trackResult(itEnds,3), 1)...
            thisDistanceToCentoid trackResult(itEnds,3) itTracks];
    end
end




saveWithHeaders([imagesFolder 'Tracks\'], 'FullInfo.dat', {'Xend', 'Yend',...
    'Xjoint','Yjoint','D2Joint','Xcentroid','Ycentroid','D2Centroid','Frame',...
    'id'}, fullEndInfo);
% save([imagesFolder 'Tracks\' fileName 'FullInfo.dat'], 'fullEndInfo', '-ASCII', '-DOUBLE');

save([imagesFolder 'Tracks\' 'FullInfo.dat'], 'fullEndInfo', '-ASCII', '-DOUBLE');
        
        

    

