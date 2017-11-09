%% Read parameters form settings file


% Author: Antoine Godin
% godin.antoine@sympatico.ca

settingsRead=inifile([imagesFolder 'Settings.ini'],'readall');

fileName     = getFromINI(settingsRead, 'fileName');
frames       = getFromINI(settingsRead, 'frames');
originPoints = getFromINI(settingsRead, 'idAndOrigin');
        % put origin and joints in the propper form
        xNear = originPoints(:,2);
        yNear = originPoints(:,3);        

%% Load data
skeletonEnds      = load([imagesFolder 'Tracks\' fileName 'Ends.dat']);
skeletonJoints    = load([imagesFolder 'Tracks\' fileName 'Joints.dat']);
fullInfo          = load([imagesFolder 'Tracks\' 'FullInfo.dat']);
trackedFrames     = unique(fullInfo(:,9));
segmentsList      = [];
phyllopodiaLength = [];

trackedFrames     = unique(fullInfo(:,9)); % column 9 is the id

warning off
itwaitbar = waitbar(0,['Making Segments by following Skeleton']);     
for itFrame=1:frames
    endsToTrack           = [];
    thisPhyllopodiaLength = [];   
    thisSkeletons = imread([imagesFolder 'nakedSkeletons\' fileName num2str(itFrame-1, '%04.f') '.tif']);
    thisFrameRows = fullInfo(find(fullInfo(:,9)==trackedFrames(itFrame)),:);
    % find origin coordinates
    counter = 1;
    for itOrigin = 1:size(originPoints, 1)
        if size(find(thisFrameRows(:,10) == originPoints(itOrigin, 1)),1)>0
            fila = thisFrameRows(find(thisFrameRows(:,10) == originPoints(itOrigin, 1)),:);
            endsToTrack(counter, 1) = fila(1,1);
            endsToTrack(counter, 2) = fila(1,2);
            endsToTrack(counter, 3) = fila(1,10);
            counter = counter+1;
        end
    end    
    rgbPath         = zeros(size(thisSkeletons,1),size(thisSkeletons,2),3);
    rgbPath(:,:,1)  = double(imread([imagesFolder fileName num2str(itFrame-1, '%04.f') '.tif']))/255;
    if size(endsToTrack,1)>0 %only does the tracing if there are ends to track in this frame
        thisFrameEnds   = [skeletonEnds(find(skeletonEnds(:,3)==itFrame),1),skeletonEnds(find(skeletonEnds(:,3)==itFrame),2)];
        thisSkeletonsWithoutJoints = thisSkeletons;
        [connectivityMatrix,russianForm,PositionVectices] = GetConnectivityPointbyPoint(thisSkeletons);
        closestPoints = getClosestJoint(xNear,yNear,PositionVectices);      
        for itPhyllopodia=1:size(endsToTrack, 1)
            idBeg = closestPoints(find(endsToTrack(itPhyllopodia,3) == originPoints(:,1)));
            idEnd = find(and(PositionVectices(:,1)==endsToTrack(itPhyllopodia,1),PositionVectices(:,2)==endsToTrack(itPhyllopodia,2)));           
            tic
            [dSP,shotestPath] = grShortPath(russianForm,idBeg,idEnd);
            toc
            shortWayImage = zeros(size(thisSkeletons));
            for itSegments=1:size(shotestPath,2)
                shortWayImage(PositionVectices(shotestPath(itSegments),2),PositionVectices(shotestPath(itSegments),1)) = 1;
            end            
            thisPhyllopodiaLength(itPhyllopodia,1) = getDistanceFollowingSkeleton(PositionVectices,shotestPath);
            thisPhyllopodiaLength(itPhyllopodia,2) = endsToTrack(itPhyllopodia,3);
            thisPhyllopodiaLength(itPhyllopodia,3) = itFrame;
            rgbPath(:,:,2) =  or(rgbPath(:,:,2),shortWayImage);
        end        
        phyllopodiaLength = [phyllopodiaLength; thisPhyllopodiaLength];
    end
    rgbPath(:,:,3) = thisSkeletons;
    imwrite(rgbPath, [imagesFolder 'FollowSkeleton\' fileName num2str(itFrame-1, '%04.f') '.jpg'], 'JPG')
    itwaitbar = waitbar(itFrame/frames,itwaitbar,['Making Segments by following Skeleton']);
end
close(itwaitbar)
saveWithHeaders([imagesFolder 'FollowSkeleton\'], 'LengthAlongPhyllopodia.dat', {'Length', 'Id','Frame'}, phyllopodiaLength);

%% Saves a log of the settings in settings@date.log
inifile([imagesFolder 'FollowSkeleton\settings@' datestr(now, 'yy.mm.dd.HH.MM') '.log'],...
    'write', settingsRead, 'tabbed');