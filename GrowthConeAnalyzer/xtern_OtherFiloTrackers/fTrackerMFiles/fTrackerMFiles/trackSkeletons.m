%% Track the ends that are the result of analyzeSkeletonMovie

% Author: Santiago Costantino 
% santiago.costantino@umontreal.ca

warning off

%% Movie Parameters

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
centroids=load([imagesFolder 'Tracks\' fileName 'CenterOfMass.dat']);

%% Track
trackResult = track(skeletonEnds, maxDisplacement, trackingParameters);
coneImage=imread([imagesFolder 'Skeletons\' fileName '0000.jpg']);

%% Make Tracks Image in a labeled-like fashion

itwaitbar = waitbar(0,'Making Tracks');
trackLinesImage=zeros(size(coneImage));
for itTrack=1:max(trackResult(:,4))
    thisTrackPixels=find(trackResult(:,4)==itTrack);
    for itLocation=1:length(thisTrackPixels)-1
        trackLinesImage=func_Drawline(trackLinesImage, ...
            trackResult(thisTrackPixels(itLocation),2),...
            trackResult(thisTrackPixels(itLocation),1),...
            trackResult(thisTrackPixels(itLocation+1),2),...
            trackResult(thisTrackPixels(itLocation+1),1), itTrack);
    end
end
%% Merge skeletons and tracks and save
itwaitbar = waitbar(1/3,itwaitbar,'Merging skeletons and tracks and save');

colores={[.5 0 0], [1 0 0], [0 .5 0], [0 1 0], [0 0 .5], [0 0 .5], ...
    [0 0 1], [1 0.5 0], [1 0 0.5], [0 1 0.5], [0.5 0.5 1], [0.5 1 0.5],...
    [1 0.5 0.5], [1 1 0.5], [1 0.5 1], [0.5 1 1], [1 1 1]};
red=zeros(size(coneImage,1), size(coneImage,2));
green=zeros(size(coneImage,1), size(coneImage,2));
blue=zeros(size(coneImage,1), size(coneImage,2));

for itTrack=1:max(trackResult(:,4))
    red(trackLinesImage==itTrack)=colores{mod(itTrack, size(colores,2))+1}(1);
    green(trackLinesImage==itTrack)=colores{mod(itTrack, size(colores,2))+1}(2);
    blue(trackLinesImage==itTrack)=colores{mod(itTrack, size(colores,2))+1}(3);
end

%% Save track movie
for itFrame=1:frames
    thisFrame=imread([imagesFolder 'Skeletons\' fileName num2str(itFrame-1, '%04.f') '.jpg']);
    redCone=thisFrame(:,:,1);
    greenCone=thisFrame(:,:,2);
    blueCone=thisFrame(:,:,3);
    redCone(red~=0)=red(red~=0)*255;
    greenCone(green~=0)=green(green~=0)*255;
    blueCone(blue~=0)=blue(blue~=0)*255;
    rbgCone(:,:,1)=redCone;
    rbgCone(:,:,2)=greenCone;
    rbgCone(:,:,3)=blueCone;
    imwrite(rbgCone, [imagesFolder '\Tracks\' fileName num2str(itFrame-1, '%04.f') '.jpg'], 'jpg')
end
    
%% Create results spreadsheet

%% Load data
skeletonEnds=load([imagesFolder 'Tracks\' fileName 'Ends.dat']);
skeletonJoints=load([imagesFolder 'Tracks\' fileName 'Joints.dat']);
skeletonCentroids=load([imagesFolder 'Tracks\' fileName 'CenterOfMass.dat']);

%% Analize tracks
itwaitbar = waitbar(2/3,itwaitbar,'Analizing tracks');

for itTracks=1:max(trackResult(:,4))
    itTracks
    for itEnds=find(trackResult(:,4)==itTracks)'
        thisEnd=[trackResult(itEnds,1), trackResult(itEnds,2)];
        theseJoints=[skeletonJoints(find(skeletonJoints(:,3)==trackResult(itEnds,3)),1),skeletonJoints(find(skeletonJoints(:,3)==trackResult(itEnds,3)),2)];
        thisTime=max(skeletonJoints(find(skeletonJoints(:,3)==2),3));
        if size(theseJoints,1) ~= 0
            [thisDistance, closestJoint]=getDistance(thisEnd, theseJoints);
            [thisDistanceToCentoid, centroid]=getDistance(thisEnd, ...
            [skeletonCentroids(thisTime, 2), skeletonCentroids(thisTime, 1)]);
            fullEndInfo=[fullEndInfo; thisEnd closestJoint thisDistance...
                    skeletonCentroids(trackResult(itEnds,3), 2) ...
                    skeletonCentroids(trackResult(itEnds,3), 1)...
            thisDistanceToCentoid trackResult(itEnds,3) itTracks];
        end
    end
end


saveWithHeaders([imagesFolder 'Tracks\'], 'FullInfo.dat', {'Xend', 'Yend',...
    'Xjoint','Yjoint','D2Joint','Xcentroid','Ycentroid','D2Centroid','Frame',...
    'id'}, fullEndInfo);
% save([imagesFolder 'Tracks\' fileName 'FullInfo.dat'], 'fullEndInfo', '-ASCII', '-DOUBLE');

save([imagesFolder 'Tracks\' 'FullInfo.dat'], 'fullEndInfo', '-ASCII', '-DOUBLE');

%% Saves a log of the settings in settings@date.log
inifile([imagesFolder 'Tracks\settings@' datestr(now, 'yy.mm.dd.HH.MM') '.log'],...
    'write', settingsRead, 'tabbed');
itwaitbar = waitbar(1,itwaitbar,'Finished with tracks');
close(itwaitbar)